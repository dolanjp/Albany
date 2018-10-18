//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//
#include <fstream>
#include "Sacado_ParameterRegistration.hpp"
#include "Albany_Utils.hpp"

#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

namespace _3DM {

//**********************************************************************
template<typename EvalT, typename Traits>
Laser_Source<EvalT, Traits>::
Laser_Source(Teuchos::ParameterList& p,
                         const Teuchos::RCP<Albany::Layouts>& dl) :
  coord_        (p.get<std::string>("Coordinate Name"),
                 dl->qp_vector),
  time        (p.get<std::string>("Time Name"),
               dl->workset_scalar),
  deltaTime   (p.get<std::string>("Delta Time Name"),
               dl->workset_scalar),
  laser_source_ (p.get<std::string>("Laser Source Name"),
                 dl->qp_scalar)
{

  this->addDependentField(coord_);
  this->addDependentField(time);
  this->addDependentField(deltaTime);
  this->addEvaluatedField(laser_source_);
 
  Teuchos::RCP<PHX::DataLayout> scalar_dl = dl->qp_scalar;
  std::vector<PHX::DataLayout::size_type> dims;
  scalar_dl->dimensions(dims);
  workset_size_ = dims[0];
  num_qps_      = dims[1];

  Teuchos::ParameterList* cond_list =
    p.get<Teuchos::ParameterList*>("Parameter List");

  Teuchos::RCP<const Teuchos::ParameterList> reflist =
    this->getValidLaser_SourceParameters();

  cond_list->validateParameters(*reflist, 0,
      Teuchos::VALIDATE_USED_ENABLED, Teuchos::VALIDATE_DEFAULTS_DISABLED);

  // dummy variable used multiple times below
  std::string type; 

    //type = cond_list->get("Laser Beam Radius Type", "Constant");
    laser_beam_radius = cond_list->get("Laser Beam Radius", 1.0);

    //type = cond_list->get("Average Laser Power Type", "Constant");
    average_laser_power = cond_list->get("Average Laser Power", 1.0);

    //type = cond_list->get("Extinction coefficient Type", "Constant");
    extinction_coeff = cond_list->get("Extinction coefficient", 1.0);

    //type = cond_list->get("Laser Pulse Frequency Type", "Constant");
    laser_pulse_frequency = cond_list->get("Laser Pulse Frequency", 1.0);

    //type = cond_list->get("Reflectivity Type", "Constant");
    reflectivity =  cond_list->get("Reflectivity", 1.0);

  this->setName("Laser_Source"+PHX::typeAsString<EvalT>());

}

//**********************************************************************
template<typename EvalT, typename Traits>
void Laser_Source<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
                      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(coord_,fm);
  this->utils.setFieldData(time,fm);
  this->utils.setFieldData(deltaTime,fm);
  this->utils.setFieldData(laser_source_,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Laser_Source<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  // current time
  const RealType t = workset.current_time;
  // time step
  ScalarT dt = deltaTime(0);
    if (dt == 0.0) dt = 1.0e-6;
  
  
  _3DM::LaserCenter Val;
  Val.t = t;
  
  RealType x, y, power_fraction;
  int power;
  LaserData_.getLaserPosition(t,Val,x,y,power,power_fraction);
  ScalarT Laser_center_x = x;
  ScalarT Laser_center_y = y;
  ScalarT Laser_power_fraction = power_fraction;


  

  // source function
  ScalarT pi = 3.1415926535897932;
  ScalarT LaserFlux_Max;
  // laser on or off
  if ( power == 1 )
    {
      LaserFlux_Max =(2.0/(pi*laser_beam_radius*laser_beam_radius))*average_laser_power*Laser_power_fraction;
    }
  else
    {
      LaserFlux_Max = 0.0;
    }
 

  

//-----------------------------------------------------------------------------------------------
  for (std::size_t cell = 0; cell < workset.numCells; ++cell) {
    for (std::size_t qp = 0; qp < num_qps_; ++qp) {
	  MeshScalarT X = coord_(cell,qp,0);
	  MeshScalarT Y = coord_(cell,qp,1);
	  MeshScalarT Z = coord_(cell,qp,2);

        ScalarT radius = sqrt((X - Laser_center_x)*(X - Laser_center_x) + (Y - Laser_center_y)*(Y - Laser_center_y));
	ScalarT depth_profile = exp(-extinction_coeff*Z);
	ScalarT surface_profile = exp(-2.0*(radius*radius)/(laser_beam_radius*laser_beam_radius));
	if (radius < laser_beam_radius)
	laser_source_(cell,qp) = extinction_coeff*(1.0 - reflectivity)*LaserFlux_Max*surface_profile*depth_profile;
	else   laser_source_(cell,qp) = 0.0;
    }
  }
}
//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<const Teuchos::ParameterList>
Laser_Source<EvalT, Traits>::
getValidLaser_SourceParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> valid_pl =
    rcp(new Teuchos::ParameterList("Valid Laser Source Params"));;
 
  //valid_pl->set<std::string>("Laser Beam Radius Type", "Constant");
  valid_pl->set<double>("Laser Beam Radius", 1.0);

  //valid_pl->set<std::string>("Average Laser Power Type", "Constant");
  valid_pl->set<double>("Average Laser Power", 1.0);

  //valid_pl->set<std::string>("Extinction coefficient Type", "Constant");
  valid_pl->set<double>("Extinction coefficient", 1.0);

  //valid_pl->set<std::string>("Laser Pulse Frequency Type", "Constant");
  valid_pl->set<double>("Laser Pulse Frequency", 1.0);

  //valid_pl->set<std::string>("Reflectivity Type", "Constant");
  valid_pl->set<double>("Reflectivity", 1.0);

  return valid_pl;
}
//**********************************************************************

}
