//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//
#include "Albany_ModelFactory.hpp"
#include "Albany_SolverFactory.hpp"
#include "MiniTensor.h"
#include "Teuchos_FancyOStream.hpp"
#include "Schwarz_Alternating.hpp"

namespace LCM {

//
//
//
SchwarzAlternating::
SchwarzAlternating(
    Teuchos::RCP<Teuchos::ParameterList> const & app_params,
    Teuchos::RCP<Teuchos::Comm<int> const> const & comm,
    Teuchos::RCP<Tpetra_Vector const> const & initial_guess)
{
  Teuchos::ParameterList &
  alt_system_params = app_params->sublist("Alternating System");

  // Get names of individual model input files
  Teuchos::Array<std::string>
  model_filenames =
      alt_system_params.get<Teuchos::Array<std::string>>("Model Input Files");

  min_iters_ = alt_system_params.get<int>("Minimum Iterations");
  max_iters_ = alt_system_params.get<int>("Maximum Iterations");
  rel_tol_ = alt_system_params.get<ST>("Relative Tolerance");
  abs_tol_ = alt_system_params.get<ST>("Absolute Tolerance");

  //number of models
  num_subdomains_ = model_filenames.size();

  // Create application name-index map used for Schwarz BC.
  Teuchos::RCP<std::map<std::string, int>>
  app_name_index_map = Teuchos::rcp(new std::map<std::string, int>);

  for (auto app_index = 0; app_index < num_subdomains_; ++app_index) {

    std::string const &
    app_name = model_filenames[app_index];

    std::pair<std::string, int>
    app_name_index = std::make_pair(app_name, app_index);

    app_name_index_map->insert(app_name_index);
  }

  //
  // Parameters
  //
  Teuchos::ParameterList &
  problem_params = app_params->sublist("Problem");

  bool const
  have_parameters = problem_params.isSublist("Parameters");

  ALBANY_ASSERT(have_parameters == false, "Parameters not supported.");

  //
  // Responses
  //
  bool const
  have_responses = problem_params.isSublist("Response Functions");

  ALBANY_ASSERT(have_responses == false, "Responses not supported.");

  //
  //
  //
  apps_.resize(num_subdomains_);
  solvers_.resize(num_subdomains_);
  convergence_ops_.resize(num_subdomains_);

  //Set up each application and model object
  for (auto m = 0; m < num_subdomains_; ++m) {

    //get parameterlist from mth model
    Albany::SolverFactory
    solver_factory(model_filenames[m], comm);

    Teuchos::ParameterList &
    params = solver_factory.getParameters();

    // Add application array for later use in Schwarz BC.
    params.set("Application Array", apps_);

    // See application index for use with Schwarz BC.
    params.set("Application Index", m);

    // App application name-index map for later use in Schwarz BC.
    params.set("Application Name Index Map", app_name_index_map);

    // Add NOX pre-post-operator for Schwarz loop convergence criterion.
    bool const
    have_piro = params.isSublist("Piro");

    ALBANY_ASSERT(have_piro == true);

    Teuchos::ParameterList &
    piro_params = params.sublist("Piro");

    bool const
    have_nox = piro_params.isSublist("NOX");

    ALBANY_ASSERT(have_nox == true);

    Teuchos::ParameterList &
    nox_params = piro_params.sublist("NOX");

    bool const
    have_solver_opts = nox_params.isSublist("Solver Options");

    ALBANY_ASSERT(have_solver_opts == true);

    Teuchos::ParameterList &
    solver_opts = nox_params.sublist("Solver Options");

    std::string const
    ppo_str{"User Defined Pre/Post Operator"};

    bool const
    have_ppo = solver_opts.isParameter(ppo_str);

    Teuchos::RCP<NOX::Abstract::PrePostOperator>
    ppo{Teuchos::null};

    if (have_ppo == true) {
      ppo = solver_opts.get<decltype(ppo)>(ppo_str);
    } else {
      ppo = Teuchos::rcp(new NOXSolverPrePostOperator);
      solver_opts.set(ppo_str, ppo);
      ALBANY_ASSERT(solver_opts.isParameter(ppo_str) == true);
    }

    Teuchos::RCP<NOXSolverPrePostOperator>
    convergence_op = Teuchos::rcp_dynamic_cast<NOXSolverPrePostOperator>(ppo);

    convergence_ops_[m] = convergence_op;

    Teuchos::RCP<Albany::Application>
    app{Teuchos::null};

    Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<ST>>
    solver = solver_factory.createAndGetAlbanyAppT(app, comm, comm);

    solvers_[m] = solver;

    apps_[m] = app;
  }

  //
  // Setup nominal values
  //
  nominal_values_ = this->createInArgsImpl();
  nominal_values_.set_x(Teuchos::null);
  nominal_values_.set_x_dot(Teuchos::null);
  nominal_values_.set_x_dot_dot(Teuchos::null);
}

//
//
//
SchwarzAlternating::
~SchwarzAlternating()
{
  return;
}

//
//
//
Teuchos::RCP<Thyra::VectorSpaceBase<ST> const>
SchwarzAlternating::
get_x_space() const
{
  return Teuchos::null;
}

//
//
//
Teuchos::RCP<Thyra::VectorSpaceBase<ST> const>
SchwarzAlternating::
get_f_space() const
{
  return Teuchos::null;
}

//
//
//
Teuchos::RCP<Thyra::VectorSpaceBase<ST> const>
SchwarzAlternating::
get_p_space(int) const
{
  return Teuchos::null;
}

//
//
//
Teuchos::RCP<Thyra::VectorSpaceBase<ST> const>
SchwarzAlternating::
get_g_space(int) const
{
  return Teuchos::null;
}

//
//
//
Teuchos::RCP<const Teuchos::Array<std::string>>
SchwarzAlternating::
get_p_names(int) const
{
  return Teuchos::null;
}

//
//
//
Teuchos::ArrayView<const std::string>
SchwarzAlternating::
get_g_names(int l) const
{
  ALBANY_ASSERT(false, "not implemented");
  return Teuchos::ArrayView<const std::string>();
}

//
//
//
Thyra::ModelEvaluatorBase::InArgs<ST>
SchwarzAlternating::
getNominalValues() const
{
  return nominal_values_;
}

//
//
//
Thyra::ModelEvaluatorBase::InArgs<ST>
SchwarzAlternating::
getLowerBounds() const
{
  return Thyra::ModelEvaluatorBase::InArgs<ST>(); // Default value
}

//
//
//
Thyra::ModelEvaluatorBase::InArgs<ST>
SchwarzAlternating::
getUpperBounds() const
{
  return Thyra::ModelEvaluatorBase::InArgs<ST>(); // Default value
}

//
//
//
Teuchos::RCP<Thyra::LinearOpBase<ST>>
SchwarzAlternating::
create_W_op() const
{
  return Teuchos::null;
}

//
//
//
Teuchos::RCP<Thyra::PreconditionerBase<ST>>
SchwarzAlternating::
create_W_prec() const
{
  return Teuchos::null;
}

//
//
//
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<ST>>
SchwarzAlternating::
get_W_factory() const
{
  return Teuchos::null;
}

//
//
//
Thyra::ModelEvaluatorBase::InArgs<ST>
SchwarzAlternating::
createInArgs() const
{
  return this->createInArgsImpl();
}

//
//
//
Teuchos::ArrayRCP<Teuchos::RCP<Albany::Application>>
SchwarzAlternating::
getApps() const
{
  return apps_;
}

//
// Create operator form of dg/dx for distributed responses
//
Teuchos::RCP<Thyra::LinearOpBase<ST>>
SchwarzAlternating::
create_DgDx_op_impl(int j) const
{
  return Teuchos::null;
}

//
// Create operator form of dg/dx_dot for distributed responses
//
Teuchos::RCP<Thyra::LinearOpBase<ST>>
SchwarzAlternating::
create_DgDx_dot_op_impl(int j) const
{
  return Teuchos::null;
}

//
// Create InArgs
//
Thyra::ModelEvaluatorBase::InArgs<ST>
SchwarzAlternating::
createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<ST>
  result;

  result.setModelEvalDescription(this->description());

  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_W_x_dot_dot_coeff, true);

  sub_inargs_.resize(num_subdomains_);
  for (auto m = 0; m < num_subdomains_; ++m) {
    sub_inargs_[m] = solvers_[m]->createInArgs();
  }

  return result;
}

//
// Create OutArgs
//
Thyra::ModelEvaluatorBase::OutArgs<ST>
SchwarzAlternating::
createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<ST>
  result;

  result.setModelEvalDescription(this->description());

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f, true);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op, true);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec, false);

  result.set_W_properties(
      Thyra::ModelEvaluatorBase::DerivativeProperties(
          Thyra::ModelEvaluatorBase::DERIV_LINEARITY_UNKNOWN,
          Thyra::ModelEvaluatorBase::DERIV_RANK_FULL,
          true));

  sub_outargs_.resize(num_subdomains_);
  for (auto m = 0; m < num_subdomains_; ++m) {
    sub_outargs_[m] = solvers_[m]->createOutArgs();
  }

  return result;
}

//
// Evaluate model on InArgs
//
void
SchwarzAlternating::
evalModelImpl(
    Thyra::ModelEvaluatorBase::InArgs<ST> const &,
    Thyra::ModelEvaluatorBase::OutArgs<ST> const &) const
{
  SchwarzLoop();
  return;
}

//
// Validate application parameters not created via a SolverFactory
// Check usage and whether necessary.
//
Teuchos::RCP<Teuchos::ParameterList const>
SchwarzAlternating::
getValidAppParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList>
  list = Teuchos::rcp(new Teuchos::ParameterList("ValidAppParams"));

  return list;
}

//
// Check usage and whether necessary.
//
Teuchos::RCP<Teuchos::ParameterList const>
SchwarzAlternating::
getValidProblemParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList>
  list = Teuchos::createParameterList("ValidProblemParams");

  return list;
}

//
// Schwarz Alternating loop
//
void
SchwarzAlternating::
SchwarzLoop() const
{
  minitensor::Vector<ST>
  norms_init(num_subdomains_, minitensor::Filler::ZEROS);

  minitensor::Vector<ST>
  norms_final(num_subdomains_, minitensor::Filler::ZEROS);

  minitensor::Vector<ST>
  norms_diff(num_subdomains_, minitensor::Filler::ZEROS);

  minitensor::Vector<int> const
  sequence(num_subdomains_, minitensor::Filler::SEQUENCE);

  int const
  iter_limit = std::max(min_iters_, max_iters_);

  std::string const
  delim(72, '=');

  Teuchos::RCP<Teuchos::FancyOStream>
  out(Teuchos::VerboseObjectBase::getDefaultOStream());

  *out << delim << '\n';
  *out << "Schwarz Alternating Method with " << num_subdomains_;
  *out << " subdomains\n";

  for (auto n = 0; n < iter_limit; ++n) {

    for (auto m = 0; m < num_subdomains_; ++m) {

      *out << delim << '\n';
      *out << "Schwarz iteration         :" << n << '\n';
      *out << "Subdomain                 :" << m << '\n';
      *out << delim << '\n';

      Thyra::ResponseOnlyModelEvaluatorBase<ST> &
      solver = *(solvers_[m]);

      Thyra::ModelEvaluatorBase::InArgs<ST> &
      in_args = sub_inargs_[m];

      Thyra::ModelEvaluatorBase::OutArgs<ST> &
      out_args = sub_outargs_[m];

      solver.evalModel(in_args, out_args);

      Teuchos::RCP<NOXSolverPrePostOperator>
      convergence_op = convergence_ops_[m];

      norms_init(m) = convergence_op->getInitialNorm();
      norms_final(m) = convergence_op->getFinalNorm();
      norms_diff(m) = convergence_op->getDifferenceNorm();
    }


    ST const
    norm_final = minitensor::norm(norms_final);

    ST const
    norm_diff = minitensor::norm(norms_diff);

    ST const
    abs_error = norm_diff;

    ST const
    rel_error = norm_final > 0.0 ? norm_diff / norm_final : norm_diff;

    *out << delim << '\n';
    *out << "Schwarz iteration         :" << n << '\n';
    *out << "Subdomain                 :" << sequence << '\n';
    *out << "Initial norms |X0|        :" << norms_init << '\n';
    *out << "Final norms   |Xf|        :" << norms_final << '\n';
    *out << "Diff norms |Xf-X0||       :" << norms_diff << '\n';
    *out << "Abs error ||Xf-X0||       :" << std::setw(24) << abs_error << '\n';
    *out << "Abs tolerance             :" << std::setw(24) << abs_tol_ << '\n';
    *out << "Rel error ||Xf-X0||/||Xf||:" << std::setw(24) << rel_error << '\n';
    *out << "Rel tolerance             :" << std::setw(24) << rel_tol_ << '\n';
    *out << delim << '\n';

    if (abs_error < abs_tol_ || rel_error < rel_tol_) {
      *out << "Schwarz loop converged.\n";
      *out << delim << '\n';
      break;
    }
  }

  return;
}

} // namespace LCM
