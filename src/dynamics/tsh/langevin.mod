	  É  ,   k820309              15.0        pWše                                                                                                           
       langevin.f90 LANGEVIN                      @                             
       #         @                                                      #RAN    #NR                                                                  
               &                                                                                                                                                    
                 
                       ð?        1.0#         @                                                     #SUB_SET_LANGEVIN_GAMMA%DABS    #GAMMA    #N_ATOM 	   #MASS 
   #G0    #DTIME                                                   DABS          D                                                     
     p          5  p        r 	       5  p        r 	                               
                                  	                    
                                  
                    
    p          5  p        r 	       5  p        r 	                               
                                       
                
  @                                    
      #         @                                                   	   #GAMMA    #MASS    #VEL_X    #VEL_Y    #VEL_Z    #GRAD_X    #GRAD_Y    #GRAD_Z    #N_ATOM             
                                                      
    p          5  p 	       r        5  p 	       r                               
                                                      
    p          5  p 	       r        5  p 	       r                               
                                                      
    p          5  p 	       r        5  p 	       r                               
                                                      
    p          5  p 	       r        5  p 	       r                               
                                                      
    p          5  p 	       r        5  p 	       r                               D                                                     
     p          5  p 	       r        5  p 	       r                               D                                                     
 	    p          5  p 	       r        5  p 	       r                               D                                                     
 
    p          5  p 	       r        5  p 	       r                                                                             #         @                                                     #SUB_LANGEVIN_NOISE%DSQRT    #GAMMA    #TEMP    #MASS    #GRAD_X    #GRAD_Y    #GRAD_Z    #N_ATOM                                                   DSQRT          
                                                      
    p          5  p        r        5  p        r                                                                      
                
                                                      
    p          5  p        r        5  p        r                               D                                                     
     p          5  p        r        5  p        r                               D                                                     
     p          5  p        r        5  p        r                               D                                                     
     p          5  p        r        5  p        r                                                                             #         @                                                        #MASS !   #VEL_X #   #VEL_Y $   #VEL_Z %   #G0 &   #TEMP '   #GRAD_X (   #GRAD_Y )   #GRAD_Z *   #N_ATOM "   #DTIME +             @                               !                    
     p          5  p 
       r "       5  p 
       r "                               @                               #                    
     p          5  p 
       r "       5  p 
       r "                               @                               $                    
     p          5  p 
       r "       5  p 
       r "                               @                               %                    
     p          5  p 
       r "       5  p 
       r "                                @                               &     
                 D @                               '     
                D @                               (                    
     p          5  p 
       r "       5  p 
       r "                              D @                               )                    
     p          5  p 
       r "       5  p 
       r "                              D @                               *                    
     p          5  p 
       r "       5  p 
       r "                               D @                               "                       @                               +     
                    fn#fn    ū   @   J   RANDOMLIB )   þ   Y       GET_STD_NORMAL+RANDOMLIB -   W     a   GET_STD_NORMAL%RAN+RANDOMLIB ,   ã  @   a   GET_STD_NORMAL%NR+RANDOMLIB    #  s       KB '            SUB_SET_LANGEVIN_GAMMA ,   3  =      SUB_SET_LANGEVIN_GAMMA%DABS -   p  ī   a   SUB_SET_LANGEVIN_GAMMA%GAMMA .   $  @   a   SUB_SET_LANGEVIN_GAMMA%N_ATOM ,   d  ī   a   SUB_SET_LANGEVIN_GAMMA%MASS *     @   a   SUB_SET_LANGEVIN_GAMMA%G0 -   X  @   a   SUB_SET_LANGEVIN_GAMMA%DTIME &     Ū       SUB_LANGEVIN_FRICTION ,   F  ī   a   SUB_LANGEVIN_FRICTION%GAMMA +   ú  ī   a   SUB_LANGEVIN_FRICTION%MASS ,   Ū  ī   a   SUB_LANGEVIN_FRICTION%VEL_X ,   b  ī   a   SUB_LANGEVIN_FRICTION%VEL_Y ,   	  ī   a   SUB_LANGEVIN_FRICTION%VEL_Z -   Ę	  ī   a   SUB_LANGEVIN_FRICTION%GRAD_X -   ~
  ī   a   SUB_LANGEVIN_FRICTION%GRAD_Y -   2  ī   a   SUB_LANGEVIN_FRICTION%GRAD_Z -   æ  @   a   SUB_LANGEVIN_FRICTION%N_ATOM #   &  ĩ       SUB_LANGEVIN_NOISE )   Û  >      SUB_LANGEVIN_NOISE%DSQRT )     ī   a   SUB_LANGEVIN_NOISE%GAMMA (   Í  @   a   SUB_LANGEVIN_NOISE%TEMP (     ī   a   SUB_LANGEVIN_NOISE%MASS *   Á  ī   a   SUB_LANGEVIN_NOISE%GRAD_X *   u  ī   a   SUB_LANGEVIN_NOISE%GRAD_Y *   )  ī   a   SUB_LANGEVIN_NOISE%GRAD_Z *   Ý  @   a   SUB_LANGEVIN_NOISE%N_ATOM $     Ā       SUB_LANGEVIN_MODIFY )   Ý  ī   a   SUB_LANGEVIN_MODIFY%MASS *     ī   a   SUB_LANGEVIN_MODIFY%VEL_X *   E  ī   a   SUB_LANGEVIN_MODIFY%VEL_Y *   ų  ī   a   SUB_LANGEVIN_MODIFY%VEL_Z '   ­  @   a   SUB_LANGEVIN_MODIFY%G0 )   í  @   a   SUB_LANGEVIN_MODIFY%TEMP +   -  ī   a   SUB_LANGEVIN_MODIFY%GRAD_X +   á  ī   a   SUB_LANGEVIN_MODIFY%GRAD_Y +     ī   a   SUB_LANGEVIN_MODIFY%GRAD_Z +   I  @   a   SUB_LANGEVIN_MODIFY%N_ATOM *     @   a   SUB_LANGEVIN_MODIFY%DTIME 