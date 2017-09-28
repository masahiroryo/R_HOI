# R_HOI
  
## **Ryo M. et al. Nonlinear higher-order abiotic interactions explain riverine biodiversity. Journal of Biogeography.**


**Corresponding Author: Masahiro Ryo (masahiroryo@gmail.com)**  
  Do not hesitate to contact me if you get in any troubles in using these scripts :)

     
         
    
### *Appendix_S2a_RF.R*   
To conduct random forest analyses with variable selection together with partial dependence plots.  
   
### *Appendix_S2b_RF_function.R*   
The function by Hapfelmeier & Ulm (2013) modified by M. Ryo for Windows OS.  
This script is called by Appendix_S2a_RF.R  
    
### *Appendix_S2c_interaction.R*   
To estimate the relative importance of interactions & to visualize 3-way interactions  
This script can run independently from Appendix_S2a_RF.R  
A dataset with 1 response and several predictors (csv file) is required.  
Ideally predictors are selected with variable selection.  
In the current version, 3-way visualization works only with 2 continuous and 1 categorial predictors  
    
### *Appendix_S2d_interaction_function.R*  
The function that was developed by M. Ryo, extending the work by Kerry & Okada (2012).  
This script is called by Appendix_S2c_interaction.R  
