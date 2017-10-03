void Compression(){

    if(Config.step==0 || fabs(rho0-Config.rhoV) > 1e-4){
        Config.step=0;
        while (rho0>Config.rhoV){

            acc = 0;   
         //   if(initCompress) exit(0);
            if(initCompress) Move_step();
            acceptance = static_cast<double>(acc)/N;
            Compression_step(); 
            std::cout << acceptance << std::endl;
        }

        std::cout <<  Config.rhoV  << " " << rho0 << std::endl;
        
        Config.rhoV = rho0;
        Config.rhoN = Config.rhoV/Config.Vsys;
        Config.Vbox = N/Config.rhoN;
        
        double ln = pow(Config.Vbox, 1.0/3.0)/Config.l[0];

        Config.l[0] *= ln;
        Config.l[1] *= ln;
        Config.l[2] *= ln;

        for( int i=0; i < Config.part.size(); i++){
            for( int vv=0; vv < 3; vv++){ 
                Config.part[i].pos[vv] *= ln;
            }
        }
        RenewList();
       
        Config.write("Save/Config.dat",1);
    }

    int zahl = 0;
    while( zahl < 10 ){ 
        acc = 0;   
        for( int k = 0; k < 10; k++) Move_step(); 
        acceptance = static_cast<double>(acc)/(10*N);
        zahl++;
        if (acceptance > 0.55){
            if( pos_lambda < maxpos) pos_lambda += 0.05;
            else pos_lambda = maxpos;
            if( ori_lambda < maxpos) ori_lambda += 0.05;
            else ori_lambda = maxpos;
            zahl=0;
        }
        if (acceptance < 0.45){
            if( pos_lambda < 0.01 ) pos_lambda = 0.01;
            else pos_lambda -= 0.01;
            if( ori_lambda < 0.01 ) ori_lambda = 0.01;
            else ori_lambda -= 0.01;
            zahl=0;
        }
        
        std::cout << acceptance << "\tpos_lambda = "<< pos_lambda << "\tand\tori_lambda = " << ori_lambda << std::endl;
    }
    Config.write("Config.dat",0);
    
}
