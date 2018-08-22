//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

// IK, 9/12/14: has Epetra_Comm! No other Epetra.

#include "Albany_ModelEvaluatorT.hpp"

#include "Albany_DistributedParameterLibrary.hpp"
#include "Albany_DistributedParameterDerivativeOpT.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Tpetra_ConfigDefs.hpp"

#include "Albany_TpetraThyraUtils.hpp"
#include "Albany_ThyraUtils.hpp"

#include "Albany_Application.hpp"

// uncomment the following to write stuff out to matrix market to debug
//#define WRITE_TO_MATRIX_MARKET

#ifdef WRITE_TO_MATRIX_MARKET
static int mm_counter_sol = 0;
static int mm_counter_res = 0;
static int mm_counter_jac = 0;
#endif  // WRITE_TO_MATRIX_MARKET

// IK, 4/24/15: adding option to write the mass matrix to matrix market file,
// which is needed
// for some applications.  Uncomment the following line to turn on.
//#define WRITE_MASS_MATRIX_TO_MM_FILE
#ifdef WRITE_MASS_MATRIX_TO_MM_FILE
#include "MatrixMarket_Tpetra.hpp"
#include "TpetraExt_MMHelpers.hpp"
#endif

Albany::ModelEvaluatorT::ModelEvaluatorT(
    const Teuchos::RCP<Albany::Application>&    app_,
    const Teuchos::RCP<Teuchos::ParameterList>& appParams)
    : app(app_),
      supports_xdot(false),
      supports_xdotdot(false),
      supplies_prec(app_->suppliesPreconditioner())
{
  Teuchos::RCP<Teuchos::FancyOStream> out =
      Teuchos::VerboseObjectBase::getDefaultOStream();

  // Parameters (e.g., for sensitivities, SG expansions, ...)
  Teuchos::ParameterList& problemParams   = appParams->sublist("Problem");
  Teuchos::ParameterList& parameterParams = problemParams.sublist("Parameters");

  num_param_vecs = parameterParams.get("Number of Parameter Vectors", 0);
  bool using_old_parameter_list = false;
  if (parameterParams.isType<int>("Number")) {
    int numParameters = parameterParams.get<int>("Number");
    if (numParameters > 0) {
      num_param_vecs           = 1;
      using_old_parameter_list = true;
    }
  }

  *out << "Number of parameter vectors  = " << num_param_vecs << std::endl;

  Teuchos::ParameterList& responseParams =
      problemParams.sublist("Response Functions");

  int  num_response_vecs       = app->getNumResponses();
  bool using_old_response_list = false;
  if (responseParams.isType<int>("Number")) {
    int numParameters = responseParams.get<int>("Number");
    if (numParameters > 0) {
      num_response_vecs       = 1;
      using_old_response_list = true;
    }
  }

  param_names.resize(num_param_vecs);
  param_lower_bds.resize(num_param_vecs);
  param_upper_bds.resize(num_param_vecs);
  for (int l = 0; l < num_param_vecs; ++l) {
    const Teuchos::ParameterList* pList =
        using_old_parameter_list ?
            &parameterParams :
            &(parameterParams.sublist(Albany::strint("Parameter Vector", l)));

    const int numParameters = pList->get<int>("Number");
    TEUCHOS_TEST_FOR_EXCEPTION(
        numParameters == 0,
        Teuchos::Exceptions::InvalidParameter,
        std::endl
            << "Error!  In Albany::ModelEvaluatorT constructor:  "
            << "Parameter vector "
            << l
            << " has zero parameters!"
            << std::endl);

    param_names[l] =
        Teuchos::rcp(new Teuchos::Array<std::string>(numParameters));
    for (int k = 0; k < numParameters; ++k) {
      (*param_names[l])[k] =
          pList->get<std::string>(Albany::strint("Parameter", k));
    }

    *out << "Number of parameters in parameter vector " << l << " = "
         << numParameters << std::endl;
  }

  Teuchos::Array<Teuchos::RCP<Teuchos::Array<std::string>>> response_names;
  response_names.resize(num_response_vecs);
  for (int l = 0; l < num_response_vecs; ++l) {
    const Teuchos::ParameterList* pList =
        using_old_response_list ?
            &responseParams :
            &(responseParams.sublist(Albany::strint("Response Vector", l)));

    bool number_exists = pList->getEntryPtr("Number");

    if (number_exists) {
      const int numParameters = pList->get<int>("Number");
      TEUCHOS_TEST_FOR_EXCEPTION(
          numParameters == 0,
          Teuchos::Exceptions::InvalidParameter,
          std::endl
              << "Error!  In Albany::ModelEvaluatorT constructor:  "
              << "Response vector "
              << l
              << " has zero parameters!"
              << std::endl);

      response_names[l] =
          Teuchos::rcp(new Teuchos::Array<std::string>(numParameters));
      for (int k = 0; k < numParameters; ++k) {
        (*response_names[l])[k] =
            pList->get<std::string>(Albany::strint("Response", k));
      }
    }
  }

  *out << "Number of response vectors  = " << num_response_vecs << std::endl;

  // Setup sacado and tpetra storage for parameters
  sacado_param_vec.resize(num_param_vecs);
  param_vecs.resize(num_param_vecs);
  param_vss.resize(num_param_vecs);
  thyra_response_vec.resize(num_response_vecs);

  Teuchos::RCP<const Teuchos::Comm<int>> commT = app->getComm();
  for (int l = 0; l < param_vecs.size(); ++l) {
    try {
      // Initialize Sacado parameter vector
      // The following call will throw, and it is often due to an incorrect
      // input line in the "Parameters" PL
      // in the input file. Give the user a hint about what might be happening
      app->getParamLib()->fillVector<PHAL::AlbanyTraits::Residual>(
          *(param_names[l]), sacado_param_vec[l]);
    } catch (const std::logic_error& le) {
      *out << "Error: exception thrown from ParamLib fillVector in file "
           << __FILE__ << " line " << __LINE__ << std::endl;
      *out << "This is probably due to something incorrect in the "
              "\"Parameters\" list in the input file, one of the lines:"
           << std::endl;
      for (int k = 0; k < param_names[l]->size(); ++k)
        *out << "      " << (*param_names[l])[k] << std::endl;

      throw le;  // rethrow to shut things down
    }

    // Create vector space for parameter vector
    Tpetra::LocalGlobal lg = Tpetra::LocallyReplicated;
    param_vss[l] = createLocallyReplicatedVectorSpace(sacado_param_vec[l].size(), commT);

    // Create Thyra vector for parameters
    param_vecs[l] = Thyra::createMember(param_vss[l]);

    Teuchos::ParameterList* pList;
    if (using_old_parameter_list) {
      pList = &parameterParams;
    } else {
      pList = &(parameterParams.sublist(Albany::strint("Parameter Vector",l)));
    }

    int numParameters = param_vss[l]->dim();

    // Loading lower bounds (if any)
    if (pList->isParameter("Lower Bounds"))
    {
      param_lower_bds[l] = Thyra::createMember(param_vss[l]);
      Teuchos::Array<ST> lb = pList->get<Teuchos::Array<ST>>("Lower Bounds");
      TEUCHOS_TEST_FOR_EXCEPTION (lb.size()!=numParameters, Teuchos::Exceptions::InvalidParameter,
                                  "Error! Input lower bounds array has the wrong dimension.\n");

      auto param_lower_bd_nonConstView = getNonconstLocalData(param_lower_bds[l]);
      for (unsigned int k = 0; k < numParameters; ++k) {
        param_lower_bd_nonConstView[k] = lb[k];
      }
    }

    // Loading upper bounds (if any)
    if (pList->isParameter("Upper Bounds"))
    {
      param_upper_bds[l] = Thyra::createMember(param_vss[l]);
      Teuchos::Array<ST> ub = pList->get<Teuchos::Array<ST>>("Upper Bounds");
      TEUCHOS_TEST_FOR_EXCEPTION (ub.size()!=numParameters, Teuchos::Exceptions::InvalidParameter,
                                  "Error! Input upper bounds array has the wrong dimension.\n");

      auto param_upper_bd_nonConstView = getNonconstLocalData(param_upper_bds[l]);
      for (unsigned int k = 0; k < numParameters; ++k) {
        param_upper_bd_nonConstView[k] = ub[k];
      }
    }

    // Loading nominal values (if any)
    auto param_vec_nonConstView = getNonconstLocalData(param_vecs[l]);
    if (pList->isParameter("Nominal Values"))
    {
      Teuchos::Array<ST> nvals = pList->get<Teuchos::Array<ST>>("Nominal Values");
      TEUCHOS_TEST_FOR_EXCEPTION (nvals.size()!=numParameters, Teuchos::Exceptions::InvalidParameter,
                                  "Error! Input nominal values array has the wrong dimension.\n");

      for (unsigned int k = 0; k < numParameters; ++k) {
        sacado_param_vec[l][k].baseValue = param_vec_nonConstView[k] = nvals[k];
      }
    } else {
      for (unsigned int k = 0; k < numParameters; ++k) {
        param_vec_nonConstView[k] = sacado_param_vec[l][k].baseValue;
      }
    }
  }

  // Setup distributed parameters
  distParamLib = app->getDistributedParameterLibrary();
  Teuchos::ParameterList& distParameterParams =
      problemParams.sublist("Distributed Parameters");
  Teuchos::ParameterList* param_list;
  num_dist_param_vecs =
      distParameterParams.get("Number of Parameter Vectors", 0);
  dist_param_names.resize(num_dist_param_vecs);
  *out << "Number of distributed parameters vectors  = " << num_dist_param_vecs
       << std::endl;
  const std::string* p_name_ptr;
  const std::string  emptyString("");
  for (int i = 0; i < num_dist_param_vecs; i++) {
    const std::string& p_sublist_name =
        Albany::strint("Distributed Parameter", i);
    param_list = distParameterParams.isSublist(p_sublist_name) ?
                     &distParameterParams.sublist(p_sublist_name) :
                     NULL;

    p_name_ptr = &distParameterParams.get<std::string>(
        Albany::strint("Parameter", i), emptyString);

    if (param_list != NULL) {
      const std::string& name_from_list =
          param_list->get<std::string>("Name", emptyString);

      p_name_ptr = (*p_name_ptr != emptyString) ? p_name_ptr : &name_from_list;

      TEUCHOS_TEST_FOR_EXCEPTION(
          (*p_name_ptr != name_from_list) && (name_from_list != emptyString),
          Teuchos::Exceptions::InvalidParameter,
          std::endl
              << "Error!  In Albany::ModelEvaluatorT constructor:  Provided "
                 "two different names for same parameter in Distributed "
                 "Parameters list: \""
              << *p_name_ptr
              << "\" and \""
              << name_from_list
              << "\""
              << std::endl);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
        !distParamLib->has(*p_name_ptr),
        Teuchos::Exceptions::InvalidParameter,
        std::endl
            << "Error!  In Albany::ModelEvaluatorT constructor:  "
            << "Invalid distributed parameter name \""
            << *p_name_ptr
            << "\""
            << std::endl);

    dist_param_names[i] = *p_name_ptr;
    // set parameters bonuds
    if (param_list) {
      Teuchos::RCP<const DistributedParameter> distParam = distParamLib->get(*p_name_ptr);
      if (param_list->isParameter("Lower Bound") &&
          (distParam->lower_bounds_vector() != Teuchos::null))
        distParam->lower_bounds_vector()->putScalar(
            param_list->get<ST>("Lower Bound"));
      if (param_list->isParameter("Upper Bound") &&
          (distParam->upper_bounds_vector() != Teuchos::null))
        distParam->upper_bounds_vector()->putScalar(
            param_list->get<ST>("Upper Bound"));
      if (param_list->isParameter("Initial Uniform Value") &&
          (distParam->vector() != Teuchos::null))
        distParam->vector()->putScalar(
            param_list->get<ST>("Initial Uniform Value"));
    }
  }

  for (int l = 0; l < app->getNumResponses(); ++l) {
    // Create Thyra vector for responses
    Teuchos::RCP<const Tpetra_Map> mapT = app->getResponse(l)->responseMapT();
    Teuchos::RCP<const Thyra_VectorSpace> gT_space =
        Thyra::createVectorSpace<ST>(mapT);
    thyra_response_vec[l] = Thyra::createMember(gT_space);
  }

  {
    // Determine the number of solution vectors (x, xdot, xdotdot)

    int num_sol_vectors =
        app->getAdaptSolMgrT()->getInitialSolution()->getNumVectors();

    if (num_sol_vectors > 1)  // have x dot
      supports_xdot = true;

    if (num_sol_vectors > 2)  // have both x dot and x dotdot
      supports_xdotdot = true;

    // Setup nominal values, lower and upper bounds
    nominalValues = this->createInArgsImpl();
    lowerBounds = this->createInArgsImpl();
    upperBounds = this->createInArgsImpl();

    // All the ME vectors are unallocated here
    allocateVectors();

    // TODO: Check if correct nominal values for parameters
    for (int l = 0; l < num_param_vecs; ++l) {
      nominalValues.set_p(l, param_vecs[l]);
      if(Teuchos::nonnull(param_lower_bds[l])) {
        lowerBounds.set_p(l, param_lower_bds[l]);
      }
      if(Teuchos::nonnull(param_upper_bds[l])) {
        upperBounds.set_p(l, param_upper_bds[l]);
      }
    }
    for (int l = 0; l < num_dist_param_vecs; ++l) {
      nominalValues.set_p(l+num_param_vecs, distParamLib->get(dist_param_names[l])->vector());
      lowerBounds.set_p(l+num_param_vecs, distParamLib->get(dist_param_names[l])->lower_bounds_vector());
      upperBounds.set_p(l+num_param_vecs, distParamLib->get(dist_param_names[l])->upper_bounds_vector());
    }
  }

  timer = Teuchos::TimeMonitor::getNewTimer("Albany: **Total Fill Time**");
}

void
Albany::ModelEvaluatorT::allocateVectors()
{
  const Teuchos::RCP<const Teuchos::Comm<int>>  commT = app->getComm();
  const Teuchos::RCP<const Thyra_VectorSpace>   vs    = app->getVectorSpace();
  const Teuchos::RCP<const Thyra_MultiVector> xMV =
      app->getAdaptSolMgrT()->getCurrentSolution();

  // Create Tpetra objects to be wrapped in Thyra

  // Create non-const versions of x_init [and x_dot_init [and x_dotdot_init]]
  const Teuchos::RCP<const Thyra_Vector> x_init = xMV->col(0);
  const Teuchos::RCP<Thyra_Vector> x_init_nonconst = x_init->clone_v();
  nominalValues.set_x(x_init_nonconst);

  // Have xdot
  if (xMV->domain()->dim() > 1) {
    const Teuchos::RCP<const Thyra_Vector> x_dot_init = xMV->col(1);
    const Teuchos::RCP<Thyra_Vector>       x_dot_init_nonconst = x_dot_init->clone_v();
    nominalValues.set_x_dot(x_dot_init_nonconst);
  }

  // Have xdotdot
  if (xMV->domain()->dim() > 2) {
    // Set xdotdot in parent class to pass to time integrator

    // GAH set x_dotdot for transient simulations. Note that xDotDot is a member
    // of Piro::TransientDecorator<ST, LO, Tpetra_GO, KokkosNode>
    const Teuchos::RCP<const Thyra_Vector> x_dotdot_init = xMV->col(2);
    const Teuchos::RCP<Thyra_Vector>       x_dotdot_init_nonconst = x_dotdot_init->clone_v();
    // IKT, 3/30/17: set x_dotdot in nominalValues for Tempus, now that
    // it is available in Thyra::ModelEvaluator
    this->xDotDot = x_dotdotT_init_nonconst;
    nominalValues.set_x_dot_dot(x_dotdot_init_nonconst);
  } else {
    this->xDotDot = Teuchos::null;
  }
}

// Overridden from Thyra::ModelEvaluator<ST>

Teuchos::RCP<const Thyra_VectorSpace>
Albany::ModelEvaluatorT::get_x_space() const
{
  Teuchos::RCP<const Tpetra_Map>                 map = app->getMapT();
  Teuchos::RCP<const Thyra_VectorSpace> x_space =
      Thyra::createVectorSpace<ST>(map);
  return x_space;
}

Teuchos::RCP<const Thyra_VectorSpace>
Albany::ModelEvaluatorT::get_f_space() const
{
  Teuchos::RCP<const Tpetra_Map>                 map = app->getMapT();
  Teuchos::RCP<const Thyra_VectorSpace> f_space =
      Thyra::createVectorSpace<ST>(map);
  return f_space;
}

Teuchos::RCP<const Thyra_VectorSpace>
Albany::ModelEvaluatorT::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_param_vecs + num_dist_param_vecs || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl
          << "Error!  Albany::ModelEvaluatorT::get_p_space():  "
          << "Invalid parameter index l = "
          << l
          << std::endl);
  Teuchos::RCP<const Tpetra_Map> map;
  if (l < num_param_vecs)
    map = tpetra_param_map[l];
  else
    map = distParamLib->get(dist_param_names[l - num_param_vecs])->map();
  Teuchos::RCP<const Thyra_VectorSpace> tpetra_param_space;
  if (map != Teuchos::null)
    tpetra_param_space = Thyra::createVectorSpace<ST>(map);
  else
    tpetra_param_space = Teuchos::null;
  return tpetra_param_space;
}

Teuchos::RCP<const Thyra_VectorSpace>
Albany::ModelEvaluatorT::get_g_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= app->getNumResponses() || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl
          << "Error!  Albany::ModelEvaluatorT::get_g_space():  "
          << "Invalid response index l = "
          << l
          << std::endl);

  return app->getResponse(l)->responseVectorSpace();
}

Teuchos::RCP<const Teuchos::Array<std::string>>
Albany::ModelEvaluatorT::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_param_vecs + num_dist_param_vecs || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl
          << "Error!  Albany::ModelEvaluatorT::get_p_names():  "
          << "Invalid parameter index l = "
          << l
          << std::endl);

  if (l < num_param_vecs) return param_names[l];
  return Teuchos::rcp(
      new Teuchos::Array<std::string>(1, dist_param_names[l - num_param_vecs]));
}

Thyra::ModelEvaluatorBase::InArgs<ST>
Albany::ModelEvaluatorT::getNominalValues() const
{
  return nominalValues;
}

Thyra::ModelEvaluatorBase::InArgs<ST>
Albany::ModelEvaluatorT::getLowerBounds() const
{
  return lowerBounds;
}

Thyra::ModelEvaluatorBase::InArgs<ST>
Albany::ModelEvaluatorT::getUpperBounds() const
{
  return upperBounds;
}

Teuchos::RCP<Thyra::LinearOpBase<ST>>
Albany::ModelEvaluatorT::create_W_op() const
{
  const Teuchos::RCP<Tpetra_Operator> W =
      Teuchos::rcp(new Tpetra_CrsMatrix(app->getJacobianGraphT()));
  return Thyra::createLinearOp(W);
}

Teuchos::RCP<Thyra::PreconditionerBase<ST>>
Albany::ModelEvaluatorT::create_W_prec() const
{
  Teuchos::RCP<Thyra::DefaultPreconditioner<ST>> W_prec =
      Teuchos::rcp(new Thyra::DefaultPreconditioner<ST>);
  Teuchos::RCP<Tpetra_Operator>         precOp = app->getPreconditionerT();
  Teuchos::RCP<Thyra::LinearOpBase<ST>> precOp_thyra =
      Thyra::createLinearOp(precOp);

  Teuchos::RCP<Thyra_LinearOp> Extra_W_op = create_W();

  W_prec->initializeRight(precOp_thyra);
  return W_prec;

  // Teko prec needs space for Jacobian as well
  // Extra_W_crs = Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(create_W(),
  // true);

  // TODO: Analog of EpetraExt::ModelEvaluator::Preconditioner does not exist in
  // Thyra yet!
  /*const bool W_prec_not_supported = true;
  TEUCHOS_TEST_FOR_EXCEPT(W_prec_not_supported);
  return Teuchos::null;*/
}

Teuchos::RCP<Thyra::LinearOpBase<ST>>
Albany::ModelEvaluatorT::create_DfDp_op_impl(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j >= num_param_vecs + num_dist_param_vecs || j < num_param_vecs,
      Teuchos::Exceptions::InvalidParameter,
      std::endl
          << "Error!  Albany::ModelEvaluatorT::create_DfDp_op_impl():  "
          << "Invalid parameter index j = "
          << j
          << std::endl);

  return Teuchos::rcp( new DistributedParameterDerivativeOpT(app, dist_param_names[j - num_param_vecs]) );
}

Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<ST>>
Albany::ModelEvaluatorT::get_W_factory() const
{
  return Teuchos::null;
}

Thyra::ModelEvaluatorBase::InArgs<ST>
Albany::ModelEvaluatorT::createInArgs() const
{
  return this->createInArgsImpl();
}

void
Albany::ModelEvaluatorT::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<ST>& finalPoint,
    const bool                                   wasSolved)
{
  // TODO
  TEUCHOS_TEST_FOR_EXCEPTION(
      true,
      Teuchos::Exceptions::InvalidParameter,
      "Calling reportFinalPoint in Albany_ModelEvaluatorT.cpp line 296"
          << std::endl);
}

Teuchos::RCP<Thyra::LinearOpBase<ST>>
Albany::ModelEvaluatorT::create_DgDx_op_impl(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j >= app->getNumResponses() || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl
          << "Error!  Albany::ModelEvaluatorT::create_DgDx_op_impl():  "
          << "Invalid response index j = "
          << j
          << std::endl);

  return Thyra::createLinearOp(app->getResponse(j)->createGradientOpT());
}

// AGS: x_dotdot time integrators not imlemented in Thyra ME yet
Teuchos::RCP<Thyra::LinearOpBase<ST>>
Albany::ModelEvaluatorT::create_DgDx_dotdot_op_impl(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j >= app->getNumResponses() || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl
          << "Error!  Albany::ModelEvaluatorT::create_DgDx_dotdot_op():  "
          << "Invalid response index j = "
          << j
          << std::endl);

  return Thyra::createLinearOp(app->getResponse(j)->createGradientOpT());
}

Teuchos::RCP<Thyra::LinearOpBase<ST>>
Albany::ModelEvaluatorT::create_DgDx_dot_op_impl(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j >= app->getNumResponses() || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl
          << "Error!  Albany::ModelEvaluatorT::create_DgDx_dot_op_impl():  "
          << "Invalid response index j = "
          << j
          << std::endl);

  return Thyra::createLinearOp(app->getResponse(j)->createGradientOpT());
}

Thyra::ModelEvaluatorBase::OutArgs<ST>
Albany::ModelEvaluatorT::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<ST> result;
  result.setModelEvalDescription(this->description());

  const int n_g = app->getNumResponses();
  result.set_Np_Ng(num_param_vecs + num_dist_param_vecs, n_g);

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f, true);

  if (supplies_prec)
    result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec, true);

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op, true);
  result.set_W_properties(Thyra::ModelEvaluatorBase::DerivativeProperties(
      Thyra::ModelEvaluatorBase::DERIV_LINEARITY_UNKNOWN,
      Thyra::ModelEvaluatorBase::DERIV_RANK_FULL,
      true));

  for (int l = 0; l < num_param_vecs; ++l) {
    result.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,
        l,
        Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL);
  }
  for (int i = 0; i < num_dist_param_vecs; i++)
    result.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,
        i + num_param_vecs,
        Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);

  for (int i = 0; i < n_g; ++i) {
    Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support;
    if (app->getResponse(i)->isScalarResponse()) {
      dgdx_support = Thyra::ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW;
    } else {
      dgdx_support = Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP;
    }

    result.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, i, dgdx_support);
    if (supports_xdot) {
      result.setSupports(
          Thyra::ModelEvaluatorBase::OUT_ARG_DgDx_dot, i, dgdx_support);
    }

    // AGS: x_dotdot time integrators not imlemented in Thyra ME yet
    // result.setSupports(
    //    Thyra::ModelEvaluatorBase::OUT_ARG_DgDx_dotdot, i, dgdx_support);

    for (int l = 0; l < num_param_vecs; l++)
      result.setSupports(
          Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
          i,
          l,
          Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL);

    if (app->getResponse(i)->isScalarResponse()) {
      for (int j = 0; j < num_dist_param_vecs; j++)
        result.setSupports(
            Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
            i,
            j + num_param_vecs,
            Thyra::ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW);
    } else {
      for (int j = 0; j < num_dist_param_vecs; j++)
        result.setSupports(
            Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
            i,
            j + num_param_vecs,
            Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
    }
  }

  return result;
}

namespace {
// As of early Jan 2015, it seems there is some conflict between Thyra's use of
// NaN to initialize certain quantities and Tpetra's v.update(alpha, x, 0)
// implementation. In the past, 0 as the third argument seemed to trigger a code
// path that does a set v <- alpha x rather than an accumulation v <- alpha x +
// beta v. Hence any(isnan(v(:))) was not a problem if beta == 0. That seems not
// to be entirely true now. For some reason, this problem occurs only in DEBUG
// builds in the sensitivities. I have not had time to fully dissect this
// problem to determine why the problem occurs only there, but the solution is
// nonetheless quite suggestive: sanitize v before calling update. I do this at
// the highest level, here, rather than in the responses.
void
sanitize_nans(const Thyra::ModelEvaluatorBase::Derivative<ST>& v)
{
  if (!v.isEmpty() && Teuchos::nonnull(v.getMultiVector()))
    ConverterT::getTpetraMultiVector(v.getMultiVector())->putScalar(0.0);
}
}  // namespace

void
Albany::ModelEvaluatorT::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<ST>&  inArgsT,
    const Thyra::ModelEvaluatorBase::OutArgs<ST>& outArgsT) const
{
#ifdef OUTPUT_TO_SCREEN
  std::cout << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  Teuchos::TimeMonitor Timer(*timer);  // start timer
  //
  // Get the input arguments
  //

  // Thyra vectors
  const Teuchos::RCP<const Thyra_Vector> x = inArgsT.get_x();
  const Teuchos::RCP<const Thyra_Vector> x_dot =
      (supports_xdot ? inArgsT.get_x_dot() : Teuchos::null);

  // Tpetra vectors
  const Teuchos::RCP<const Tpetra_Vector> xT = Albany::getConstTpetraVector(x);
  const Teuchos::RCP<const Tpetra_Vector> x_dotT = Albany::getConstTpetraVector(x_dot);

  // IKT, 3/30/17: the following logic is meant to support both the Thyra
  // time-integrators in Piro
  //(e.g., trapezoidal rule) and the second order time-integrators in Tempus.
  Teuchos::RCP<Tpetra_Vector> x_dotdotT;
  ST                          omega;
  if (supports_xdotdot == true) {
    omega = inArgsT.get_W_x_dot_dot_coeff();
    // The following case is to support second order time-integrators in Piro
    if (std::abs(omega) < 1.0e-14) {
      if (Teuchos::nonnull(this->get_x_dotdot())) {
        x_dotdotT = ConverterT::getTpetraVector(this->get_x_dotdot());
        omega     = this->get_omega();
      } else {
        x_dotdotT = Teuchos::null;
        omega     = 0.0;
      }
    }
    // The following case is for second-order time-integrators in Tempus
    else {
      if (inArgsT.supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot)) {
        Teuchos::RCP<const Tpetra_Vector> x_dotdotT_temp =
            Albany::getConstTpetraVector(inArgsT.get_x_dot_dot());
        x_dotdotT = Teuchos::rcp(new Tpetra_Vector(*x_dotdotT_temp));
      } else {
        x_dotdotT = Teuchos::null;
        omega     = 0.0;
      }
    }
  } else {
    x_dotdotT = Teuchos::null;
    omega     = 0.0;
  }
  const Teuchos::RCP<const Thyra_Vector> x_dotdot = Albany::createConstThyraVector(x_dotdotT);

  const ST alpha = (Teuchos::nonnull(x_dotT) || Teuchos::nonnull(x_dotdotT)) ?
                       inArgsT.get_alpha() :
                       0.0;
  const ST beta = (Teuchos::nonnull(x_dotT) || Teuchos::nonnull(x_dotdotT)) ?
                      inArgsT.get_beta() :
                      1.0;

  bool const is_dynamic =
      Teuchos::nonnull(x_dotT) || Teuchos::nonnull(x_dotdotT);

#if defined(ALBANY_LCM)
  ST const curr_time = is_dynamic == true ? inArgsT.get_t() : getCurrentTime();
#else
  const ST curr_time =
      (Teuchos::nonnull(x_dotT) || Teuchos::nonnull(x_dotdotT)) ?
          inArgsT.get_t() :
          0.0;
#endif  // ALBANY_LCM

  double dt = 0.0; //time step 
  if (is_dynamic == true) {
    dt = inArgsT.get_step_size(); 
  }

  for (int l = 0; l < num_param_vecs+num_dist_param_vecs; ++l) {
    const Teuchos::RCP<const Thyra_Vector> p = inArgsT.get_p(l);
    if (Teuchos::nonnull(p)) {

      if(l<num_param_vecs){
        auto p_constView = getLocalData(p);
        ParamVec& sacado_param_vector = sacado_param_vec[l];
        for (unsigned int k = 0; k < sacado_param_vector.size(); ++k)
          sacado_param_vector[k].baseValue = p_constView[k];
      } else {
        distParamLib->get(dist_param_names[l-num_num_param_vecs])->vector()->assign(p);
      }
    }
  }

  //
  // Get the output arguments
  //
  auto f_out = outArgsT.get_f();
  auto W_op_out = outArgsT.get_W_op();

  /*
   * Commenting this out for now, cause it seems it is never used. If that turns out to
   * be the case, then remove altogether.
   *
   * // Get preconditioner operator, if requested
   * Teuchos::RCP<Tpetra_Operator> WPrec_out;
   * if (outArgsT.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
   *   // IKT, 12/19/16: need to verify that this is correct
   *   // LB, 7/25/18: this is currently pointless because: a) we're hiding WPrec_out
   *   //              from outside the if with a temporary, and b) we never set the
   *   //              WPrec_out in the outArgs.
   *   Teuchos::RCP<Tpetra_Operator> WPrec_out = app->getPreconditionerT();
   * }
   */

#ifdef WRITE_STIFFNESS_MATRIX_TO_MM_FILE
    // IK, 4/24/15: write stiffness matrix to matrix market file
    // Warning: to read this in to MATLAB correctly, code must be run in serial.
    // Otherwise Mass will have a distributed Map which would also need to be
    // read in to MATLAB for proper reading in of Mass.
    // IMPORTANT NOTE: keep this call BEFORE the computation of the actual jacobian,
    //                 so you don't overwrite the jacobian.
    app->computeGlobalJacobian(
        0.0, 1.0, 0.0, curr_time,
        x, x_dot, x_dotdot,
        sacado_param_vec,
        Teuchos::null, W_op_out);

    writeMatrixMarket(W_op_out,"stiffness.mm");
    writeMatrixMarket(W_op_out->range(),"range_space.mm");
    writeMatrixMarket(W_op_out->domain(),"domain_space.mm");
#endif

#ifdef WRITE_MASS_MATRIX_TO_MM_FILE
    // IK, 4/24/15: write mass matrix to matrix market file
    // Warning: to read this in to MATLAB correctly, code must be run in serial.
    // Otherwise Mass will have a distributed Map which would also need to be
    // read in to MATLAB for proper reading in of Mass.
    // IMPORTANT NOTE: keep this call BEFORE the computation of the actual jacobian,
    //                 so you don't overwrite the jacobian.
    app->computeGlobalJacobian(
        1.0, 0.0, 0.0, curr_time,
        x, x_dot, x_dotdot,
        sacado_param_vec,
        Teuchos::null, W_op_out);

    writeMatrixMarket(W_op_out,"mass.mm");
    writeMatrixMarket(W_op_out->range(),"range_space.mm");
    writeMatrixMarket(W_op_out->domain(),"domain_space.mm");
#endif

  //
  // Compute the functions
  //
  bool f_already_computed = false;

  // W matrix
  if (Teuchos::nonnull(W_op_out)) {
    app->computeGlobalJacobian(
        alpha, beta, omega, curr_time,
        x, x_dot, x_dotdot,
        sacado_param_vec,
        f_out, W_op_out, dt);
    f_already_computed = true;
  }
  /*
   *  Commenting this out for now, cause it seems it is never used. If that turns out to
   *  be the case, then remove altogether.
   *
   *  if (Teuchos::nonnull(WPrec_out)) {
   *    app->computeGlobalJacobian(
   *        alpha, beta, omega, curr_time,
   *        x, x_dot, x_dotdot,
   *        sacado_param_vec,
   *        f_out, Extra_W_op);
   *    f_already_computed = true;
   *
   *    app->computeGlobalPreconditionerT(Extra_W_crs, WPrec_out);
   *  }
   */

  // df/dp
  for (int l = 0; l < num_param_vecs; ++l) {
    const Teuchos::RCP<Thyra::MultiVectorBase<ST>> dfdp_out =
        outArgsT.get_DfDp(l).getMultiVector();

    if (Teuchos::nonnull(dfdp_out)) {
      const Teuchos::RCP<ParamVec> p_vec =
          Teuchos::rcpFromRef(sacado_param_vec[l]);

      app->computeGlobalTangent(
          0.0, 0.0, 0.0, curr_time, false,
          x, x_dot, x_dotdot, sacado_param_vec, p_vec.get(),
          Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null,
          f_out, Teuchos::null, dfdp_out);

      f_already_computed = true;
    }
  }

  // distributed df/dp
  for (int i=0; i<num_dist_param_vecs; i++) {
	  const Teuchos::RCP<Thyra_LinearOp> dfdp_out = outArgsT.get_DfDp(i+num_param_vecs).getLinearOp();
    if (dfdp_out != Teuchos::null) {
      Teuchos::RCP<DistributedParameterDerivativeOpT> dfdp_op =
        Teuchos::rcp_dynamic_cast<DistributedParameterDerivativeOpT>(dfdp_out);
      dfdp_op->set(curr_time, x, x_dot, x_dotdot,
                   Teuchos::rcp(&sacado_param_vec,false));
    }
  }

  // f
  if (app->is_adjoint) {
    const Thyra::ModelEvaluatorBase::Derivative<ST> f_derivT(
        outArgsT.get_f(), Thyra::ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW);

    const Thyra::ModelEvaluatorBase::Derivative<ST> dummy_derivT;

    const int response_index = 0;  // need to add capability for sending this in
    app->evaluateResponseDerivative(
        response_index, curr_time,
        x, x_dot, x_dotdot,
        sacado_param_vec,
        NULL,
        Teuchos::null,
        f_derivT,
        dummy_derivT,
        dummy_derivT,
        dummy_derivT);
  } else if (Teuchos::nonnull(f_out) && !f_already_computed) {
    app->computeGlobalResidual(
        curr_time,
        x, x_dot, x_dotdot,
        sacado_param_vec,
        f_out,
        dt);
  }

  // Response functions
  for (int j = 0; j < outArgsT.Ng(); ++j) {
    Teuchos::RCP<Thyra_Vector> g_out = outArgsT.get_g(j);

    const Thyra::ModelEvaluatorBase::Derivative<ST> dgdxT_out = outArgsT.get_DgDx(j);
    Thyra::ModelEvaluatorBase::Derivative<ST> dgdxdotT_out;

    if (supports_xdot) {
      dgdxdotT_out = outArgsT.get_DgDx_dot(j);
    }

    //    const Thyra::ModelEvaluatorBase::Derivative<ST> dgdxdotdotT_out =
    //    this->get_DgDx_dotdot(j);
    const Thyra::ModelEvaluatorBase::Derivative<ST> dgdxdotdotT_out;

    sanitize_nans(dgdxT_out);
    sanitize_nans(dgdxdotT_out);
    sanitize_nans(dgdxdotdotT_out);

    // dg/dx, dg/dxdot
    if (!dgdxT_out.isEmpty() || !dgdxdotT_out.isEmpty()) {
      const Thyra::ModelEvaluatorBase::Derivative<ST> dummy_derivT;
      app->evaluateResponseDerivative(
          j, curr_time, x, x_dot, x_dotdot,
          sacado_param_vec,
          NULL,
          g_out,
          dgdxT_out,
          dgdxdotT_out,
          dgdxdotdotT_out,
          dummy_derivT);
      // Set g_out to null to indicate that g_out was evaluated.
      g_out = Teuchos::null;
    }

    // dg/dp
    for (int l = 0; l < num_param_vecs; ++l) {
      const Teuchos::RCP<Thyra::MultiVectorBase<ST>> dgdp_out =
          outArgsT.get_DgDp(j, l).getMultiVector();

      if (Teuchos::nonnull(dgdp_out)) {
        const Teuchos::RCP<ParamVec> p_vec = Teuchos::rcpFromRef(sacado_param_vec[l]);

        app->evaluateResponseTangent(
            j, alpha, beta, omega, curr_time, false,
            x, x_dot, x_dotdot, sacado_param_vec, p_vec.get(),
            Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null,
            g_out, Teuchos::null, dgdp_out);
        // Set g_out to null to indicate that g_out was evaluated.
        g_out = Teuchos::null;
      }
    }

    // Need to handle dg/dp for distributed p
    for (int l = 0; l < num_dist_param_vecs; l++) {
      const Teuchos::RCP<Thyra::MultiVectorBase<ST>> dgdp_out =
          outArgsT.get_DgDp(j, l + num_param_vecs).getMultiVector();

      if (Teuchos::nonnull(dgdp_out)) {
        dgdp_out->assign(0.);
        app->evaluateResponseDistParamDeriv(
            j, curr_time, x, x_dot, x_dotdot,
            sacado_param_vec,
            dist_param_names[l],
            dgdp_out);
      }
    }

    if (Teuchos::nonnull(g_out)) {
      app->evaluateResponse(j, curr_time, x, x_dot, x_dotdot, sacado_param_vec, g_out);
    }
  }

#ifdef WRITE_TO_MATRIX_MARKET
  Albany::writeMatrixMarket(x, "sol", mm_counter_sol);
  ++mm_counter_sol;
#endif

#ifdef WRITE_TO_MATRIX_MARKET
  Albany::writeMatrixMarket(f_out, "res", mm_counter_res);
  ++mm_counter_res;
#endif

#ifdef WRITE_TO_MATRIX_MARKET
  Albany::writeMatrixMarket(W_op_out, "jac", mm_counter_jac);
  ++mm_counter_jac;
#endif
}

Thyra::ModelEvaluatorBase::InArgs<ST>
Albany::ModelEvaluatorT::createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<ST> result;
  result.setModelEvalDescription(this->description());

  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x, true);

#if defined(ALBANY_LCM)
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_step_size, true);
#endif

  if (supports_xdot) {
    result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot, true);
    result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t, true);
    result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_step_size, true);
    result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha, true);
    result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta, true);
  }

  if (supports_xdotdot) {
    result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot, true);
    result.setSupports(
        Thyra::ModelEvaluatorBase::IN_ARG_W_x_dot_dot_coeff, true);
  }
  result.set_Np(num_param_vecs + num_dist_param_vecs);

  return result;
}
