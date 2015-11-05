3drimc : 3DRISM MC code

install:

   - download
   - setup path:
         export PATH=$PATH:$(pwd)
         echo 'export PATH=$PATH:'$(pwd) >> ~/.bashrc

   - compile 3drimc_realpath: 
          
        gcc -o 3drimc_realpath 3drimc_realpath.c

   - create "hosts"
        mkdir hosts
        mkdir hosts/proc1
        ...
        mkdir hosts/procN

    - run the simulation

       3drimc_run_on_host hosts/proc1 states rismmol_files 100 10

     - you can run the last command in parallel (changing hosts/proc1 to hosts/proc2, 3 etc )
     - you can qsub that command
       

    

