void RenewList(){ // Initialise Neighbour List

    Config.W[0] = (int)(Config.l[0]/rlist); 
    Config.W[1] = (int)(Config.l[1]/rlist); 
    Config.W[2] = (int)(Config.l[2]/rlist); 
    
    Config.w[0] = Config.W[0]/Config.l[0];
    Config.w[1] = Config.W[1]/Config.l[1];
    Config.w[2] = Config.W[2]/Config.l[2];
	
    Config.head.resize(Config.W[0]*Config.W[1]*Config.W[2]); 
	Config.link.resize(Config.part.size()); 

	for( int i=0; i<Config.head.size(); i++) Config.head[i] = -1;
	for( int i=0; i<Config.link.size(); i++) Config.link[i] = -1;

	for( int i=0; i<Config.part.size(); i++){ 
	    
        sx = Config.part[i].pos[0]*Config.w[0];
        if(sx==Config.W[0]) sx = Config.W[0] - 1; 
        sy = Config.part[i].pos[1]*Config.w[1]; 
        if(sy==Config.W[1]) sy = Config.W[1] - 1; 
        sz = Config.part[i].pos[2]*Config.w[2]; 
        if(sz==Config.W[2]) sz = Config.W[2] - 1; 
   
        int k = sx + Config.W[0]*(sy + Config.W[1]*sz); 

        Config.part[i].cell = k;
        if(Config.head[k]==-1) Config.head[k] = i;
        else{
           Config.link[i] = Config.head[k];
           Config.head[k] = i; 
        }
    }

};
