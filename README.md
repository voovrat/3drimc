3drimc : 3DRISM MC code

Requirements:
  - python
  - octave
  - RISM-MOL-Tools ( github.com/voovrat/RISM-MOL-Tools )
  - voovrat-utils ( github.com/voovrat/voovrat-utils ) 

install:

   - download
   - setup paths:
     
         export PATH=$PATH:$(pwd)

         echo 'export PATH=$PATH:'$(pwd) >> ~/.bashrc
         
         echo "path(path,'"$(pwd)"');" >> ~/.octaverc


   - compile 3drimc_realpath: 
          
        gcc -o 3drimc_realpath 3drimc_realpath.c

   - create "hosts"
   - 
        mkdir hosts
        mkdir hosts/proc1
        ...
        mkdir hosts/procN

    - run the simulation

       3drimc_run_on_host hosts/proc1 states rismmol_files 100 10

     - you can run the last command in parallel (changing hosts/proc1 to hosts/proc2, 3 etc )
     - you can qsub that command
       

    

