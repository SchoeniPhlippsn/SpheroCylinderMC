bool overlapSPH ( class teilchen sphere1, class teilchen sphere2, std::vector<double> l){

    std::vector<double> R = DeltaR(sphere1.pos,sphere2.pos,l);
	double Rsq = scal_p(R,R);

    
    if(Rsq > (sphere1.r+sphere2.r)*(sphere1.r+sphere2.r) ) return false;
    else return true;
}

bool overlapSPSPH ( class teilchen cone, class teilchen sphere, std::vector<double> l){
	
    std::vector<double> R = DeltaR(sphere.pos,cone.pos,l);
    double Rsq = scal_p(R,R);
    if(Rsq > (sphere.r+cone.Rdist+cone.R)*(sphere.r+cone.Rdist+cone.R) ) return false;
    else{
        double Rw = scal_p(R,cone.ori);
        double Rw2 = Rw*Rw;
        double norm = Rsq - Rw2;
        if(norm > (sphere.r+cone.R)*(sphere.r+cone.R) ) return false;
        else{
            if(norm < 1e-6 ) return true; 
            else{
                norm = 1/sqrt(norm);

                std::vector<double> wT=R;

                for (int v = 0; v < 3; v++){
                    wT[v] -= Rw*cone.ori[v];
                    wT[v] *= norm;
                }

                double RwT = scal_p(R,wT);

                double RR = Rsq - 2*RwT*cone.a + cone.a*cone.a;
                double B = Rw - RwT*cone.c + cone.a*cone.c;
                double D = 1 + cone.c*cone.c;
            
                double mu;
                double distance = sphere.r;
                mu = B/D;

                if(mu>cone.rdist || mu<-cone.Rdist){
                    if(Rw > 0 ){
                        mu = cone.rdist;
                        distance += cone.r;
                    }else{
                        mu = -cone.Rdist;
                        distance += cone.R;
                    }
                    RR = Rsq - 2*mu*Rw + mu*mu;
                    if(RR < distance*distance) return true;
                    else return false;
                }else{
                    RR = RR - 2*mu*B + mu*mu*D;
                    if(RR < distance*distance) return true;
                    else return false;
                }
            }
        }
    }
}

bool overlapSP ( class teilchen cone1, class teilchen cone2, std::vector<double> l){
	
    std::vector<double> R = DeltaR(cone2.pos,cone1.pos,l);
    double Rsq = scal_p(R,R);
    if(Rsq > 4*(cone1.Rdist+cone1.R)*(cone1.Rdist+cone1.R)) return false;
    else{
        class teilchen sphere=cone1;
        sphere.r = cone1.R;        
        for (int v = 0; v < 3; v++) sphere.pos[v] -= cone1.Rdist*cone1.ori[v];
        if(overlapSPSPH(cone2,sphere,l)) return true;
        else{
            sphere=cone1;
            for (int v = 0; v < 3; v++) sphere.pos[v] += cone1.rdist*cone1.ori[v];
            if(overlapSPSPH(cone2,sphere,l)) return true;
            else{
                sphere=cone2;
                sphere.r = cone2.R;        
                for (int v = 0; v < 3; v++) sphere.pos[v] -= cone2.Rdist*cone2.ori[v];
                if(overlapSPSPH(cone1,sphere,l)) return true;
                else{
                    sphere=cone2;
                    for (int v = 0; v < 3; v++) sphere.pos[v] += cone2.rdist*cone2.ori[v];
                    if(overlapSPSPH(cone1,sphere,l)) return true;
                    else{
                        double Rw1 = scal_p(R,cone1.ori);
                        double Rw2 = scal_p(R,cone2.ori);

                        std::vector<double> R_neu = R;
                        double w1w2 = scal_p(cone1.ori,cone2.ori);

                        double lambdaa = 0;
                        double muu = 0;

                        if( fabs(w1w2*w1w2 - 1) < 1e-6 )  return false;
                        else{
                            if(fabs(w1w2) < 1e-6){
                                lambdaa = -Rw2;
                                muu = Rw1;
                                if(lambdaa > cone2.rdist) lambdaa = cone2.rdist;
                                else if(muu < -cone1.Rdist) muu = -cone1.Rdist;
                                if(muu > cone1.rdist) muu = cone1.rdist;
                                else if(lambdaa < -cone2.Rdist) lambdaa = -cone2.Rdist;
                            }else{
                                lambdaa = (Rw2 - Rw1*w1w2)/(w1w2*w1w2-1);
                                if(lambdaa > cone2.rdist) lambdaa = cone2.rdist;
                                if(lambdaa < -cone2.Rdist) lambdaa = -cone2.Rdist;
                                muu = (Rw2 + lambdaa)/w1w2;
                                if(muu > cone1.rdist){ 
                                    muu = cone1.rdist;
                                    lambdaa= muu*w1w2-Rw2;
                                    if(lambdaa > cone2.rdist) lambdaa = cone2.rdist;
                                    else if(lambdaa < -cone2.Rdist) lambdaa = -cone2.Rdist;
                                }else{
                                    if(muu < -cone1.Rdist){ 
                                        muu = -cone1.Rdist;
                                        lambdaa= muu*w1w2-Rw2;
                                        if(lambdaa > cone2.rdist) lambdaa = cone2.rdist;
                                        else if(lambdaa < -cone2.Rdist) lambdaa = -cone2.Rdist;
                                    }
                                }
                            }
                            for (int v = 0; v < 3; v++) R_neu[v] = R[v] + lambdaa*cone2.ori[v] - muu*cone1.ori[v];

                            double Rsq_neu = scal_p(R_neu,R_neu);
                            double R_neuw1 = scal_p(R_neu,cone1.ori);
                            double R_neuw2 = scal_p(R_neu,cone2.ori);

                            double norm1 = Rsq_neu-R_neuw1*R_neuw1;
                            double norm2 = Rsq_neu-R_neuw2*R_neuw2;
                            
                            if( norm1 < 1e-6 || norm2 < 1e-6 ) return false;
                            else{
                                norm1 = 1/sqrt(norm1);
                                norm2 = 1/sqrt(norm2);

                                std::vector<double> w1T=R_neu;
                                std::vector<double> w2T=R_neu;


                                for (int v = 0; v < 3; v++){
                                    w1T[v] -= R_neuw1*cone1.ori[v];
                                    w1T[v] *= norm1;
                                    w2T[v] -= R_neuw2*cone2.ori[v];
                                    w2T[v] *= -norm2;
                                }


                                std::vector<double> w1T_1 = kreuz(cone1.ori,w1T,cone2.ori);
                                std::vector<double> w2T_1 = kreuz(cone2.ori,w2T,cone1.ori);

                                double normal1=0;
                                double normal2=0;

                                for (int v = 0; v < 3; v++){
                                    w1T[v] += 0.25*(w1w2+1)*w1T_1[v];
                                    normal1 += w1T[v]*w1T[v];
                                    w2T[v] += 0.25*(w1w2+1)*w2T_1[v];
                                    normal2 += w2T[v]*w2T[v];
                                }

                                normal1 = 1/sqrt(normal1);
                                normal2 = 1/sqrt(normal2);

                                for (int v = 0; v < 3; v++){
                                    w1T[v] *= normal1;
                                    w2T[v] *= normal2;
                                }

                                
                                double Rw1T = scal_p(R,w1T);
                                double Rw2T = scal_p(R,w2T);
                                double w1Tw2T = scal_p(w1T,w2T);
                                double w1w2T = scal_p(cone1.ori,w2T);
                                double w1Tw2 = scal_p(w1T,cone2.ori);

                                double A = Rw2 - Rw2T*cone2.c - cone2.c*cone2.a - cone2.a*w1Tw2 + cone2.a*cone2.c*w1Tw2T;
                                double B = -Rw1 + Rw1T*cone2.c - cone2.a*cone2.c - cone2.a*w1w2T + cone2.a*cone2.c*w1Tw2T;
                                double C = 1 + cone2.c*cone2.c;
                                double D = 1 + cone2.c*cone2.c;
                                double E = -w1w2 + w1Tw2*cone2.c + w1w2T*cone2.c - w1Tw2T*cone2.c*cone2.c;

                                
                                double lambda, mu;

                                if( fabs(D*C - E*E) < 1e-6 ){
                                    if(Rw2 < 0){ 
                                        lambda = cone2.rdist;        
                                    }else{ 
                                        lambda = -cone2.Rdist;
                                    }

                                    mu = -(B+lambda*E)/D;

                                    if(mu>cone1.rdist || mu<-cone1.Rdist) return false;
                                    else{ 
                                        double RR = Rsq + 2*Rw2T*cone2.a - 2*Rw1T*cone2.a + cone2.a*cone2.a - 2*w1Tw2T*cone2.a*cone2.a + cone2.a*cone2.a;
                                        RR = RR + 2*lambda*A + 2*mu*B + lambda*lambda*C + mu*mu*D + 2*mu*lambda*E;
                                        if(RR < 1e-6 ) return true;
                                        else return false; 
                                    }
                                }else{
                                    lambda = (E*B-A*D)/(D*C-E*E);
                                    if(lambda>cone2.rdist || lambda<-cone2.Rdist) return false;
                                    else{
                                        mu = -(B+lambda*E)/D;

                                        if(mu>cone1.rdist || mu<-cone1.Rdist) return false;
                                        else{
                                            double RR = Rsq + 2*Rw2T*cone2.a - 2*Rw1T*cone2.a + cone2.a*cone2.a - 2*w1Tw2T*cone2.a*cone2.a + cone2.a*cone2.a;
                                            RR = RR + 2*lambda*A + 2*mu*B + lambda*lambda*C + mu*mu*D + 2*mu*lambda*E;

                                            double antipara = Rw1T + lambda*w1Tw2 - lambda*cone1.c*w1Tw2T + cone1.a*w1Tw2T - cone1.a + mu*cone1.c;

                                            if(RR < 1e-8 || antipara < 0 ) return true;
                                            else return false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
