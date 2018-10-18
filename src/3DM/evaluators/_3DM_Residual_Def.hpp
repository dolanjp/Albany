//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"
#include "PHAL_Utilities.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

namespace _3DM {

  //**********************************************************************
  template<typename EvalT, typename Traits>
  Phase_Residual<EvalT, Traits>::
  Phase_Residual(const Teuchos::ParameterList& p,
		 const Teuchos::RCP<Albany::Layouts>& dl) :
    w_bf_             (p.get<std::string>("Weighted BF Name"), dl->node_qp_scalar),
    w_grad_bf_        (p.get<std::string>("Weighted Gradient BF Name"), dl->node_qp_vector),
    T_                (p.get<std::string>("Temperature Name"), dl->qp_scalar),
    T_grad_           (p.get<std::string>("Temperature Gradient Name"), dl->qp_vector),
    //T_dot_            (p.get<std::string>("Temperature Time Derivative Name"), dl->qp_scalar),
    k_                (p.get<std::string>("Thermal Conductivity Name"), dl->qp_scalar),
    rho_cp_           (p.get<std::string>("rho_Cp Name"), dl->qp_scalar),
    source_           (p.get<std::string>("Source Name"), dl->qp_scalar),
    laser_source_     (p.get<std::string>("Laser Source Name"), dl->qp_scalar),
    time              (p.get<std::string>("Time Name"), dl->workset_scalar),
    psi1_             (p.get<std::string>("Psi1 Name"), dl->qp_scalar),
    psi2_             (p.get<std::string>("Psi2 Name"), dl->qp_scalar),
    phi1_             (p.get<std::string>("Phi1 Name"), dl->qp_scalar),
    phi2_             (p.get<std::string>("Phi2 Name"), dl->qp_scalar),
    energyDot_        (p.get<std::string>("Energy Rate Name"), dl->qp_scalar),
    deltaTime         (p.get<std::string>("Delta Time Name"), dl->workset_scalar),
    residual_         (p.get<std::string>("Residual Name"), dl->node_scalar)
  {
    this->addDependentField(w_bf_);
    this->addDependentField(w_grad_bf_);
    this->addDependentField(T_);
    this->addDependentField(T_grad_);
    //this->addDependentField(T_dot_);
    this->addDependentField(k_);
    this->addDependentField(rho_cp_);
    this->addDependentField(source_);
    this->addDependentField(laser_source_);
    this->addDependentField(phi1_);
    this->addDependentField(phi2_);
    this->addDependentField(psi1_);
    this->addDependentField(psi2_);
    this->addDependentField(energyDot_);
    this->addDependentField(time);
    this->addDependentField(deltaTime);

    this->addEvaluatedField(residual_);
  
  
  

  
    std::vector<PHX::Device::size_type> dims;
    w_grad_bf_.fieldTag().dataLayout().dimensions(dims);
    workset_size_ = dims[0];
    num_nodes_    = dims[1];
    num_qps_      = dims[2];
    num_dims_     = dims[3];

    Temperature_Name_ = p.get<std::string>("Temperature Name")+"_old";


    this->setName("Phase_Residual"+PHX::typeAsString<EvalT>());
  }

  //**********************************************************************
  template<typename EvalT, typename Traits>
  void Phase_Residual<EvalT, Traits>::
  postRegistrationSetup(typename Traits::SetupData d,
			PHX::FieldManager<Traits>& fm)
  {
    this->utils.setFieldData(w_bf_,fm);
    this->utils.setFieldData(w_grad_bf_,fm);
    this->utils.setFieldData(T_,fm);
    this->utils.setFieldData(T_grad_,fm);
    //this->utils.setFieldData(T_dot_,fm);
    this->utils.setFieldData(k_,fm);
    this->utils.setFieldData(rho_cp_,fm);
    this->utils.setFieldData(source_,fm);
    this->utils.setFieldData(laser_source_,fm);
    this->utils.setFieldData(time,fm);
    this->utils.setFieldData(deltaTime,fm);
    this->utils.setFieldData(phi1_,fm);
    this->utils.setFieldData(phi2_,fm);
    this->utils.setFieldData(psi1_,fm);
    this->utils.setFieldData(psi2_,fm);
    this->utils.setFieldData(energyDot_,fm);
    this->utils.setFieldData(residual_,fm);


    term1_ = Kokkos::createDynRankView(k_.get_view(), "term1_", workset_size_,num_qps_,num_dims_);
  }

  //**********************************************************************
  template<typename EvalT, typename Traits>
  void Phase_Residual<EvalT, Traits>::
  evaluateFields(typename Traits::EvalData workset)
  {
    // time step
    ScalarT dt = deltaTime(0);
    typedef Intrepid2::FunctionSpaceTools<PHX::Device> FST;

  //  if (dt == 0.0) dt = 1.0e-15;
    //grab old temperature
    Albany::MDArray T_old = (*workset.stateArrayPtr)[Temperature_Name_];
    
 /*   // Compute Temp rate
    for (std::size_t cell = 0; cell < workset.numCells; ++cell)
      {
        for (std::size_t qp = 0; qp < num_qps_; ++qp)
	  {
            T_dot_(cell, qp) = (T_(cell, qp) - T_old(cell, qp)) / dt;
	  }
      }
 */
    // diffusive term
    FST::scalarMultiplyDataData<ScalarT> (term1_, k_.get_view(), T_grad_.get_view());
    
    // zero out residual
    for (int cell = 0; cell < workset.numCells; ++cell) {
      for (int node = 0; node < num_nodes_; ++node) {
        residual_(cell,node) = 0.0;
      }
    }


    //THESE ARE HARD CODED NOW. NEEDS TO BE CHANGED TO USER INPUT LATER
    ScalarT Coeff_volExp = 65.2e-6; //per kelvins
    ScalarT Ini_temp = 300; //kelvins
   
 
      for (int cell = 0; cell < workset.numCells; ++cell) {
	for (int qp = 0; qp < num_qps_; ++qp) {
	  for (int node = 0; node < num_nodes_; ++node) {
	    residual_(cell, node) += (
				      w_grad_bf_(cell, node, qp, 0) * term1_(cell, qp, 0)
				      + w_grad_bf_(cell, node, qp, 1) * term1_(cell, qp, 1)
				      + w_grad_bf_(cell, node, qp, 2) * term1_(cell, qp, 2));
	  }
	}
      }
      // heat source from laser 
      for (int cell = 0; cell < workset.numCells; ++cell) {
	for (int qp = 0; qp < num_qps_; ++qp) {
	  for (int node = 0; node < num_nodes_; ++node) {
	    residual_(cell, node) -= (w_bf_(cell, node, qp) * laser_source_(cell, qp));
	  }
	}
      }
      // all other problem sources
      for (int cell = 0; cell < workset.numCells; ++cell) {
	for (int qp = 0; qp < num_qps_; ++qp) {
	  for (int node = 0; node < num_nodes_; ++node) {
	    residual_(cell, node) -= (w_bf_(cell, node, qp) * source_(cell, qp));
	  }
	}
      }
      // transient term
      for (int cell = 0; cell < workset.numCells; ++cell) {
	for (int qp = 0; qp < num_qps_; ++qp) {
	  for (int node = 0; node < num_nodes_; ++node) {
	    residual_(cell, node) += (w_bf_(cell, node, qp) * energyDot_(cell, qp));
	  }
	}
      }
	  
    }
         
  

  //*********************************************************************
}
