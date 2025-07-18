# Specify solver headers
set( solidMechanicsSolvers_headers
     ${solidMechanicsSolvers_headers}
     contact/ContactSolverBase.hpp
     contact/ContactFields.hpp
     contact/LogLevelsInfo.hpp
     contact/SolidMechanicsEmbeddedFractures.hpp
     contact/SolidMechanicsLagrangeContact.hpp
     contact/SolidMechanicsLagrangeContactBubbleStab.hpp
     contact/SolidMechanicsAugmentedLagrangianContact.hpp
     contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp
     contact/kernels/SolidMechanicsDisplacementJumpUpdateKernels.hpp
     contact/kernels/SolidMechanicsEFEMKernelsBase.hpp
     contact/kernels/SolidMechanicsEFEMKernels.hpp
     contact/kernels/SolidMechanicsEFEMStaticCondensationKernels.hpp
     contact/kernels/SolidMechanicsEFEMKernelsHelper.hpp
     contact/kernels/SolidMechanicsALMKernelsBase.hpp
     contact/kernels/SolidMechanicsALMKernels.hpp
     contact/kernels/SolidMechanicsConformingContactKernelsHelper.hpp
     contact/kernels/SolidMechanicsContactFaceBubbleKernels.hpp
     contact/kernels/SolidMechanicsLagrangeContactKernels.hpp 
     PARENT_SCOPE )


# Specify solver sources
set( solidMechanicsSolvers_sources
     ${solidMechanicsSolvers_sources}
     contact/ContactSolverBase.cpp
     contact/SolidMechanicsEmbeddedFractures.cpp
     contact/SolidMechanicsLagrangeContact.cpp
     contact/SolidMechanicsLagrangeContactBubbleStab.cpp
     contact/SolidMechanicsAugmentedLagrangianContact.cpp
     PARENT_SCOPE )