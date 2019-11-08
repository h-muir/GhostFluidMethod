//Basic Ghost Fluid method implemented with HLLC solver for the Euler equations
//commented nicely for assignments

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <assert.h>
using namespace std;

typedef vector<double> Vector;

template<typename T>
void printVector(vector<T> vec){
    int len = vec.size();
    for(int i; i<len; ++i){
        cout << vec[i] << endl;
    }
}//vector-array print function templated to print for any vector type

//three class definitions for: 
//   - Primitive variables (W) 
//   - Conservative variables (U)
//   - Conservative fluxes (F) 

class Prim{
    public:
        Prim(double, double, double, double);
        double rho;
        double u;
        double p;
        double gamma;
};
Prim::Prim(double _rho, double _u, double _p, double _gamma){
    rho = _rho; u = _u; p=_p; gamma=_gamma;
}

class ConsU{
    public:
        ConsU(double, double, double, double);
        double rho;
        double rhou;
        double E;
        double gamma;
};
ConsU::ConsU(double _rho, double _rhou, double _E, double _gamma){
    rho = _rho; rhou = _rhou; E = _E; gamma = _gamma;
}

class ConsF{
    public:
        ConsF(double, double, double);
        double mass;
        double mom;
        double en;
};
ConsF::ConsF(double _mass, double _mom, double _en){
    mass = _mass; mom = _mom; en = _en;
}


//Print functions for classes-----------------------------------------------:
void printW(Prim W){
    printf("---------------\n W: \n rho: %f \n u:   %f \n p:   %f \n e:   %f \n gamma: %f \n --------------\n",\
            W.rho, W.u, W.p, W.p/(W.rho*(W.gamma-1)), W.gamma);
}

void printU(ConsU U){
    printf("---------------\n U: \n rho:  %f \n rhou: %f \n E:    %f \n gamma:    %f \n --------------\n"\
            ,U.rho, U.rhou, U.E, U.gamma);
}

void printF(ConsF F){
    printf("---------------\n F: \n mass:  %f \n mom: %f \n en:    %f \n --------------\n"\
            ,F.mass, F.mom, F.en);
}
//--------------------------------------------------------------------------/


//Conversion functions:

ConsU primtoconsU(Prim W){
    double rho = W.rho;
    double rhou = W.rho*W.u;
    double E = 0.5*W.rho*pow(W.u,2) + W.p/(W.gamma-1);
    double gamma = W.gamma;
    ConsU U(rho, rhou, E, gamma);
    return U;
}

ConsF primtoconsF(Prim W){
    double mass = W.rho*W.u;
    double mom = W.rho*pow(W.u,2) + W.p;
    double E = 0.5*W.rho*pow(W.u,2) + W.p/(W.gamma-1);
    double en = W.u*(E + W.p);
    ConsF F(mass, mom, en);
    return F;
}

Prim consUtoprim(ConsU U){
    double rho = U.rho;
    double u = U.rhou/rho;
    double p = (U.gamma-1)*(U.E - 0.5*rho*pow(u,2));
    double gamma = U.gamma;
    Prim W(rho, u, p, gamma);
    return W;
}


#include "TestCases_GhostFluid.H"

ConsU delta(double omega, ConsU U_l, ConsU U_m, ConsU U_n, string scheme){
    /* omega = [-1, 1]
     * U_l = U_{i-1}, U_m = Ui, U_n = U{i+1}
     */
    ConsU delta_i_l(U_m.rho - U_l.rho, U_m.rhou - U_l.rhou, U_m.E - U_l.E, 0);
    ConsU delta_i_n(U_n.rho - U_m.rho, U_n.rhou - U_m.rhou, U_n.E - U_m.E, 0);
    
    //alternative conventions for dividing by zero RHS slope:
    // -3 encoded to represent -ve infinity
    // +3 encoded to represent +ve infinity 
    double r_rho = delta_i_l.rho;
    if(delta_i_n.rho == 0){
        if(delta_i_l.rho == 0){
            r_rho = 0;           //both left and right slopes are zero
        }else if(delta_i_l.rho < 0){
            r_rho = -3;          //negative slope on left, zero slope on right
        }else{
            r_rho = 3;           //positive slope on left, zero slope on right
        }
    }else{r_rho = r_rho/delta_i_n.rho;}     //permissible division by delta_i_n
    double r_rhou =  delta_i_l.rhou;
    if(delta_i_n.rhou == 0){
       if(delta_i_l.rhou == 0){
            r_rhou = 0;          //both left and right slopes are zero
        }else if(delta_i_l.rhou < 0){
            r_rhou = -3;         //negative slope on left, zero slope on right
        }else{
            r_rhou = 3;          //positive slope on left, zero slope on right
        }
    }else{r_rhou = r_rhou/delta_i_n.rhou;}  //permissible division by delta_i_n
    double r_E = delta_i_l.E;
    if(delta_i_n.E == 0){ 
        if(delta_i_l.E == 0){
            r_E = 0;          //both left and right slopes are zero
        }else if(delta_i_l.E < 0){
            r_E = -3;         //negative slope on left, zero slope on right
        }else{
            r_E = 3;          //positive slope on left, zero slope on right
        }
    }else{r_E = r_E/delta_i_n.E;}           //permissible division by delta_i_n

    ConsU delta_i(0.5*(1+omega)*delta_i_l.rho+0.5*(1-omega)*delta_i_n.rho,
                  0.5*(1+omega)*delta_i_l.rhou+0.5*(1-omega)*delta_i_n.rhou,
                  0.5*(1+omega)*delta_i_l.E+0.5*(1-omega)*delta_i_n.E, 0);

    ConsU sigma_R(2/(1-omega+(1+omega)*r_rho), 2/(1-omega+(1+omega)*r_rhou),
                  2/(1-omega+(1+omega)*r_E), 0);
    
    //adjustments for infinity encoded cases:
    if(r_rho == -3 || r_rho == 3){
        sigma_R.rho = 0;
    }
    if(r_rhou == -3 || r_rhou == 3){
        sigma_R.rhou = 0;
    }
    if(r_E == -3 || r_E == 3){
        sigma_R.E = 0;
    }

    if(scheme == "SuperBee"){
        ConsU delta_bar(0,0,0,0);
        if(r_rho < 0){
            delta_bar.rho = 0;
        }else if(0 <= r_rho && r_rho < 0.5){
            delta_bar.rho = 2*r_rho*delta_i.rho;
        }else if(0.5 <= r_rho && r_rho < 1){
            delta_bar.rho = 1*delta_i.rho;
        }else if(1 <= r_rho && r_rho < 2){
            delta_bar.rho = ((r_rho < sigma_R.rho)? 
                    r_rho*delta_i.rho : sigma_R.rho*delta_i.rho);
        }else{
            delta_bar.rho = ((sigma_R.rho < 2)? 
                    sigma_R.rho*delta_i.rho : 2*delta_i.rho);
;
        }
        if(r_rhou < 0){
            delta_bar.rhou = 0;
        }else if(0 <= r_rhou && r_rhou < 0.5){
            delta_bar.rhou = 2*r_rhou*delta_i.rhou;
        }else if(0.5 <= r_rhou && r_rhou < 1){
            delta_bar.rhou = 1*delta_i.rhou;
        }else if(1 <= r_rhou && r_rhou < 2){
            delta_bar.rhou = ((r_rhou < sigma_R.rhou)? 
                    r_rhou*delta_i.rhou : sigma_R.rhou*delta_i.rhou);
        }else{
            delta_bar.rhou = ((sigma_R.rhou < 2)? 
                    sigma_R.rhou*delta_i.rhou : 2*delta_i.rhou);
        }
        if(r_E < 0){
            delta_bar.E = 0;
        }else if(0 <= r_E && r_E < 0.5){
            delta_bar.E = 2*r_E*delta_i.E;
        }else if(0.5 <= r_E && r_E < 1){
            delta_bar.E = 1*delta_i.E;
        }else if(1 <= r_E && r_E < 2){
            delta_bar.E = ((r_E < sigma_R.E)? 
                    r_E*delta_i.E : sigma_R.E*delta_i.E);
        }else{
            delta_bar.E = ((sigma_R.E < 2)? 
                    sigma_R.E*delta_i.E : 2*delta_i.E);
        }

        return delta_bar;

    }else{

        return delta_i;

    }
}


class Cell{
    public:
        double aL, aR;                  //left and right sound speeds
        double SL, SR, Sstar, Splus;    //wave speeds
        ConsU UL, UR, ULstar, URstar;   //4 cons. states
        ConsF FL, FR, FLstar, FRstar;   //4 cons. fluxes
        ConsU U;                        //tbd state
        ConsF Fout;                     //flux out (F_i+1/2)
        
        //must declare constructor here:
        Cell(Prim WL, Prim WR):
            //must initialize all class objects here:
            UL(0,0,0,0), UR(0,0,0,0), ULstar(0,0,0,0), URstar(0,0,0,0),
            FL(0,0,0), FR(0,0,0), FLstar(0,0,0), FRstar(0,0,0),
            U(0,0,0,0), Fout(0,0,0) {
            //direct wave speed estimates:
            aL = pow(fabs(WL.gamma*WL.p/WL.rho), 0.5);
            aR = pow(fabs(WR.gamma*WR.p/WR.rho), 0.5);
            SL = ((WL.u - aL < WR.u - aR)? WL.u - aL : WR.u - aR);
            SR = ((WL.u + aL > WR.u + aR)? WL.u + aL : WR.u + aR);
            Sstar = (WR.p - WL.p + WL.rho*WL.u*(SL-WL.u) - WR.rho*WR.u*(SR-WR.u))/\
                    (WL.rho*(SL-WL.u)-WR.rho*(SR-WR.u));
            Splus = ((fabs(WL.u)+aL > fabs(WR.u)+aR)? fabs(WL.u)+aL : fabs(WR.u)+aR);
            //W-->U and W-->F conversions:
            UL = primtoconsU(WL);
            UR = primtoconsU(WR);
            FL = primtoconsF(WL);
            FR = primtoconsF(WR);
            //Star region calculations:
            ULstar.rho = WL.rho*((SL-WL.u)/(SL-Sstar));
            ULstar.rhou = ULstar.rho*Sstar;
            ULstar.E = ULstar.rho*(UL.E/WL.rho+(Sstar-WL.u)*(Sstar + WL.p/(WL.rho*(SL-WL.u))));
            ULstar.gamma = WL.gamma;
            //assert(ULstar.E >=0); //check that total energy is non-negative
            URstar.rho = WR.rho*((SR-WR.u)/(SR-Sstar));
            URstar.rhou = URstar.rho*Sstar;
            URstar.E = URstar.rho*(UR.E/WR.rho+(Sstar-WR.u)*(Sstar + WR.p/(WR.rho*(SR-WR.u))));
            URstar.gamma = WR.gamma;
            //assert(URstar.E >=0); //check that total energy is non-negative
            FLstar.mass = FL.mass + SL*(ULstar.rho - UL.rho);
            FLstar.mom = FL.mom + SL*(ULstar.rhou - UL.rhou);
            FLstar.en = FL.en + SL*(ULstar.E - UL.E);
            FRstar.mass = FR.mass + SR*(URstar.rho - UR.rho);
            FRstar.mom = FR.mom + SR*(URstar.rhou - UR.rhou);
            FRstar.en = FR.en + SR*(URstar.E - UR.E);
            }

};

vector<Prim> updateW(vector<Prim> W_vec, int n, double &dt, double dx, double CFL, double &Smax,
        bool MUSCL, string scheme){

    //MUSCL-Hancock parameters:
    ConsU UL(0,0,0,0);            
    ConsU Ul(0,0,0,0);
    ConsU Ui(0,0,0,0);
    ConsU Un(0,0,0,0);
    ConsU UR(0,0,0,0);
    ConsU UbarL(0,0,0,0);
    ConsU UbarR(0,0,0,0);
    Prim WbarL(0,0,0,0);
    Prim WbarR(0,0,0,0);
    ConsF FbarL(0,0,0);
    ConsF FbarR(0,0,0);    
    double omega = 0;
    ConsU delta_i(0,0,0,0);
    vector<Prim> WbarL_vec(n, W_vec[0]);
    vector<Prim> WbarR_vec(n, W_vec[0]);

    //update paramters 
    vector<Prim> W_vec_new = W_vec; //initialised variable for subsequent updating  
    Cell cell0(W_vec[0], W_vec[0]);
    vector<Cell> cell_vec(n, cell0);
    ConsU Unew(0,0,0,0);
    Prim Wnew(0,0,0,0);

    for(int i = 0; i<n; ++i){
        //MUSCL-Hancock steps 1+2:
        if(i == 0){
            Ul = primtoconsU(W_vec[i]);
        }else{
            Ul = primtoconsU(W_vec[i-1]); //left-side boundary
        }

        Ui = primtoconsU(W_vec[i]);

        if(i == n-1){
            Un = primtoconsU(W_vec[i]);
        }else{
            Un = primtoconsU(W_vec[i+1]); //right-side boundary
        }
        
        
        if(MUSCL == 1){        
            delta_i = delta(omega, Ul, Ui, Un, scheme);
        }

        UL = ConsU(Ui.rho - 0.5*delta_i.rho,
                   Ui.rhou - 0.5*delta_i.rhou,
                   Ui.E - 0.5*delta_i.E,
                   Ui.gamma);
        UR = ConsU(Ui.rho + 0.5*delta_i.rho,
                   Ui.rhou + 0.5*delta_i.rhou,
                   Ui.E + 0.5*delta_i.E,
                   Ui.gamma);
        WbarL = consUtoprim(UL);
        WbarR = consUtoprim(UR);
        FbarL = primtoconsF(WbarL);
        FbarR = primtoconsF(WbarR);
        UbarL.rho = UL.rho + 0.5*(dt/dx)*(FbarL.mass - FbarR.mass);
        UbarL.rhou = UL.rhou + 0.5*(dt/dx)*(FbarL.mom - FbarR.mom);
        UbarL.E = UL.E + 0.5*(dt/dx)*(FbarL.en - FbarR.en);
        UbarL.gamma = UL.gamma;
        UbarR.rho = UR.rho + 0.5*(dt/dx)*(FbarL.mass - FbarR.mass);
        UbarR.rhou = UR.rhou + 0.5*(dt/dx)*(FbarL.mom - FbarR.mom);
        UbarR.E = UR.E + 0.5*(dt/dx)*(FbarL.en - FbarR.en);
        UbarR.gamma = UR.gamma;

        WbarL = consUtoprim(UbarL);
        WbarR = consUtoprim(UbarR);
        
        WbarL_vec[i] = WbarL;
        WbarR_vec[i] = WbarR;
        

    }
    for(int i = 0; i<n; ++i){
        if(i==n-1){
            cell_vec[i] = Cell(WbarR_vec[i], WbarR_vec[i]);
        }else{
            cell_vec[i] = Cell(WbarR_vec[i], WbarL_vec[i+1]);                
        }

        if(cell_vec[i].Splus > Smax){
            Smax = cell_vec[i].Splus;
        }
    }
    //some printing options for progress updates during long simulations: 
    //cout << "nT: " << nT << endl;
    //cout << "dt: " << dt << endl;
    //cout << "Smax: " << Smax << endl;
    assert(Smax > 0);
    dt = CFL*dx/Smax;   //dt update
    assert(dt > 0);     //check this is true rather than enforcing fabs(dt)
    Cell cell = cell0;  //any initiation
    Cell prev_cell = cell0;
    Cell zerocell = cell0;
    Cell onecell = cell0;

    for(int i = 0; i<n; ++i){

        cell = cell_vec[i];
        Ui = primtoconsU(W_vec[i]);

        if( 0 < cell.SL){
            cell.U = cell.UL;
            cell.Fout = cell.FL;
        }else if(cell.SL <= 0 && 0 < cell.Sstar){
            cell.U = cell.ULstar;
            cell.Fout = cell.FLstar;
        }else if(cell.Sstar <= 0 && 0 < cell.SR){
            cell.U = cell.URstar;
            cell.Fout = cell.FRstar;
        }else if(cell.SR <= 0){
            cell.U = cell.UR;
            cell.Fout = cell.FR;
        }
        
        if(i == 0){
            prev_cell = cell;
            zerocell = cell;
        }
        if(i == 1){
            onecell = cell;
        }
        
        //update
        Unew.rho = Ui.rho + dt/dx * (prev_cell.Fout.mass - cell.Fout.mass);
        Unew.rhou = Ui.rhou + dt/dx * (prev_cell.Fout.mom - cell.Fout.mom);
        Unew.E = Ui.E + dt/dx * (prev_cell.Fout.en - cell.Fout.en);
        Unew.gamma = Ui.gamma;
        Wnew = consUtoprim(Unew);
        W_vec_new[i] = Wnew;
        prev_cell = cell; 
                    
    }
    
    //LEFT side boundary update due to RHS defined flux: F(U(xi+1/2))
    cell = cell_vec[0];
    Ui = primtoconsU(W_vec[0]);
    Unew.rho = Ui.rho + dt/dx * (zerocell.Fout.mass - onecell.Fout.mass);
    Unew.rhou = Ui.rhou + dt/dx * (zerocell.Fout.mom - onecell.Fout.mom);
    Unew.E = Ui.E + dt/dx * (zerocell.Fout.en - onecell.Fout.en);
    Unew.gamma = Ui.gamma;
    Wnew = consUtoprim(Unew);
    W_vec_new[0] = Wnew;

    return W_vec_new;
}

Vector update_level_set(Vector level_set, vector<Prim> W1_vec, vector<Prim> W2_vec, 
        int n, double dt, double dx, bool re_init){
    Vector updated_level_set(n);
    double Dxf;   //forward difference;
    double Dxb;   //backward difference
    int index_c;  //index of last cell before contact wave
    double eps_x; //subcell distance for linear interpolation 
    for(int i = 0; i<n; ++i){
        if(i == 0 || i == n-1){
            if(level_set[i] <= 0){
                updated_level_set[i] = level_set[i] - W1_vec[i].u*dt;
            }else{
                updated_level_set[i] = level_set[i] - W2_vec[i].u*dt;
            }
        }else{
            Dxf = level_set[i+1] - level_set[i];
            Dxb = level_set[i] - level_set[i-1];
            if(level_set[i] <= 0){
                index_c = i;
                if(W1_vec[i].u <= 0){
                    updated_level_set[i] = level_set[i] - W1_vec[i].u*(dt/dx)*Dxf;
                }else{
                    updated_level_set[i] = level_set[i] - W1_vec[i].u*(dt/dx)*Dxb;
                }
            }else{
                if(W2_vec[i].u <= 0){
                    updated_level_set[i] = level_set[i] - W2_vec[i].u*(dt/dx)*Dxf;
                }else{
                    updated_level_set[i] = level_set[i] - W2_vec[i].u*(dt/dx)*Dxb;
                }
            }
        }
    }
    //linear interpolation to exact location of contact wave:
    eps_x = updated_level_set[index_c]*dx/
            (updated_level_set[index_c] - updated_level_set[index_c+1]);

    if(re_init == 1){ //linear re-initialisation
        for(int i = 0; i < n; ++i){
            updated_level_set[i] = i*dx - (index_c*dx + eps_x); //signed distance to contact wave
        }
    }
    return updated_level_set;
}

void Casper(vector<Prim> &W1_vec, vector<Prim> &W2_vec, int index_LG, int n, int nG,
        bool isobar_fix, bool set_all_G){
    //Casper updates our u, p and rho(entropy condition) variables
    //stealthy update by reference using pointers: &
    //index_LG is the index of the last real left fluid cell
    //nG is the chosen number of ghost cells
    assert((index_LG - nG) > 0);
    assert((index_LG + nG) < n);
    double rhoG1, rhoG2;          //for entropy condition calculated rho
    int index_RG = index_LG + 1;  //the first real fluid cell in material 2 (right of ghost cells)

    if(isobar_fix == 1){
        rhoG1 = pow(((W1_vec[index_LG].p * pow(W1_vec[index_LG-1].rho,W1_vec[index_LG-1].gamma))/
                W1_vec[index_LG-1].p), 1/W1_vec[index_LG-1].gamma);
        W1_vec[index_LG].rho = rhoG1;

        rhoG2 = pow(((W2_vec[index_RG].p * pow(W2_vec[index_RG+1].rho,W2_vec[index_RG+1].gamma))/
                W2_vec[index_RG+1].p), 1/W2_vec[index_RG+1].gamma);
        W2_vec[index_RG].rho = rhoG2; 
    }  

    int left_marker;
    int right_marker;

    if(set_all_G == 0){
        left_marker = index_LG - (nG-1);
        right_marker = index_RG + nG;
    }else if(set_all_G == 1){
        left_marker = 0;
        right_marker = n;
    }

    for(int i = left_marker; i <= index_LG; ++i){
        W2_vec[i].u = W1_vec[i].u;
        W2_vec[i].p = W1_vec[i].p;
        rhoG2 = pow(((W2_vec[i].p * pow(W2_vec[index_RG].rho,W2_vec[index_RG].gamma))/
                W2_vec[index_RG].p), 1/W2_vec[index_RG].gamma); 
        //formula: [Pg*rhoI^gammaI/PI]^1/gammaI
        W2_vec[i].rho = rhoG2;
    }
    for(int i=index_RG; i < right_marker; ++i){
        W1_vec[i].u = W2_vec[i].u;
        W1_vec[i].p = W2_vec[i].p;
        rhoG1 = pow(((W1_vec[i].p * pow(W1_vec[index_LG].rho,W1_vec[index_LG].gamma))/
                W1_vec[index_LG].p), 1/W1_vec[index_LG].gamma); 
        W1_vec[i].rho = rhoG1;
   }
}


vector< vector<Prim> > generate_solution(TestCase CASE, Vector time_vec, int nG, bool set_all_G,
                                            bool isobar_fix, bool MUSCL, string scheme, 
                                            bool re_init){
    int n = CASE.n;
    double CFL = CASE.CFL;
    double T0 = CASE.T0;
    double Tf = CASE.Tf;
    double dx = CASE.dx;
    double T = T0;
    double L = CASE.L;
    double x0 = CASE.x0;
    double x1 = CASE.x1;
    double x2 = CASE.x2;
    Prim W1_0 = CASE.W1;
    Prim W2_0 = CASE.W2;
    Prim W3_0 = CASE.W3;
    Prim W4_0 = CASE.W4;
    int num_states = CASE.num_states;
    Vector space = CASE.space;
    Vector level_set = CASE.level_set;
    Vector level_set2 = CASE.level_set2;
    Vector level_set3 = CASE.level_set3;
    Vector updated_level_set(n); 
    Vector updated_level_set2(n); 
    Vector updated_level_set3(n); 

    Cell cell0(W1_0, W1_0);
    vector<Cell> cell_vec1(n, cell0);
    vector<Cell> cell_vec2(n, cell0);

    double Smax=0;
    double dt = 0;

    vector< vector<Prim> > MATRIX1;     //solution is generated as a space*time matrix 
    vector< vector<Prim> > MATRIX2;     //of the primitive variables

    vector< vector<Prim> > MATRIXreal;

    vector<Prim> W1_vec(n, W1_0);  
    vector<Prim> W2_vec(n, W2_0);   
    vector<Prim> W1_vec_new = W1_vec;   //initialised variable for subsequent updating 
    vector<Prim> W2_vec_new = W2_vec;   //initialised variable for subsequent updating 
    vector<Prim> W_vec_real(n, W1_0);
    vector<Prim> W_vec_ghost(n, W1_0);
    
    //additional initialised vectors for 3 or 4 states cases
    vector<Cell> cell_vec3(n, cell0);
    vector< vector<Prim> > MATRIX3;
    vector<Prim> W3_vec(n, W3_0);
    vector<Prim> W3_vec_new = W3_vec;
    vector<Cell> cell_vec4(n, cell0);
    vector< vector<Prim> > MATRIX4;
    vector<Prim> W4_vec(n, W4_0);
    vector<Prim> W4_vec_new = W4_vec;

    int index_LG;                       //Index of the last real fluid cell on the left 
                                        //(LG: left of ghost)
    int index_LG2;                      //for states 2-3 interface
    int index_LG3;                      //for states 3-4 interface
    
    //zero time initialisation:
    for(int i=0; i < n; ++i){ 
        if(level_set[i] <= 0){
            W_vec_real[i] = W1_vec[i];
            index_LG = i;
        }else{
            W_vec_real[i] = W2_vec[i];
        }
    }
    if(num_states >= 3){
        assert(x1 > x0);
        for(int i = 0; i<n; ++i){
            if(level_set2[i] <= 0){
                index_LG2 = i;
            }else{
                W_vec_real[i] = W3_vec[i];
            }
        }
    }
    if(num_states >= 4){
        assert(x2 > x1);
        for(int i = 0; i<n; ++i){
            if(level_set3[i] <= 0){
                index_LG3 = i;
            }else{
                W_vec_real[i] = W4_vec[i];
            }
        }
    }

    MATRIXreal.push_back(W_vec_real);
    MATRIX1.push_back(W1_vec);
    MATRIX2.push_back(W2_vec);
    if(num_states >= 3){MATRIX3.push_back(W3_vec);}
    if(num_states >= 4){MATRIX4.push_back(W4_vec);}

    for(int i = 1; i < n; ++i){
        cell0 = Cell(W1_vec[i-1], W1_vec[i]);
        if(cell0.Splus > Smax){
            Smax = cell0.Splus;
        }
        cell0 = Cell(W2_vec[i-1], W2_vec[i]);
        if(cell0.Splus > Smax){
            Smax = cell0.Splus;
        }
        if(num_states >= 3){
            cell0 = Cell(W3_vec[i-1], W3_vec[i]);
            if(cell0.Splus > Smax){
                Smax = cell0.Splus;
            }
        }
        if(num_states >= 4){
            cell0 = Cell(W4_vec[i-1], W4_vec[i]);
            if(cell0.Splus > Smax){
                Smax = cell0.Splus;
            }
        }

    }
    //global Smax calculated to take smallest of the stable timesteps between 
    //the two or three materials
    assert(Smax > 0);
    dt = CFL*dx/Smax;
    double dt1;
    double dt2;
    double dt3;
    time_vec.push_back(0);
    int nT = 1; // time step counter - zero time is the first time step
                    
    while(T < Tf){
        //Casper fixes the ghost cell properties
        Casper(W1_vec, W2_vec, index_LG, n, nG, isobar_fix, set_all_G);
        if(num_states >= 3){
            Casper(W2_vec, W3_vec, index_LG2, n, nG, isobar_fix, set_all_G); 
        }
        if(num_states >= 4){
            Casper(W3_vec, W4_vec, index_LG3, n, nG, isobar_fix, set_all_G); 
        }             
        Smax = 0;
        W1_vec_new = updateW(W1_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
        dt1 = dt;
        W2_vec_new = updateW(W2_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
        dt2 = dt;
        if(dt < dt1){
            //first time step was larger than second: 
            //re-do first material update with smaller time-step
            W1_vec_new = updateW(W1_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
        }
        if(num_states >= 3){
            W3_vec_new = updateW(W3_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
            dt3 = dt;
            if(dt < dt2){
                W1_vec_new = updateW(W1_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
                W2_vec_new = updateW(W2_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
                //recalculate with smallest timestep of the three states
                //&dt is passed by reference and only updated if Smax increases 
                //(ie. dt decreases)
            }
        }
        if(num_states >= 4){
            W4_vec_new = updateW(W4_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
            if(dt < dt3){
                W1_vec_new = updateW(W1_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
                W2_vec_new = updateW(W2_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
                W3_vec_new = updateW(W3_vec, n, dt, dx, CFL, Smax, MUSCL, scheme);
                //recalculate with smallest timestep of the four states
            }
        } 
        updated_level_set = update_level_set(level_set, W1_vec_new, W2_vec_new, 
                n, dt, dx, re_init);
        if(num_states >= 3){
            updated_level_set2 = update_level_set(level_set2, W2_vec_new, W3_vec_new, 
                    n, dt, dx, re_init);
        }
        if(num_states >= 4){
            updated_level_set3 = update_level_set(level_set3, W3_vec_new, W4_vec_new, 
                    n, dt, dx, re_init);
        }
        for(int i=0; i < n; ++i){ 
            if(updated_level_set[i] <= 0){
                W_vec_real[i] = W1_vec_new[i];
                index_LG = i;
            }else{
                W_vec_real[i] = W2_vec_new[i];
            }
        }
        if(num_states >= 3){
            for(int i = 0; i<n; ++i){
                if(updated_level_set2[i] <= 0){
                    index_LG2 = i;
                }else{
                    W_vec_real[i] = W3_vec[i];
                }
            }
        }
        if(num_states >= 4){
            for(int i = 0; i<n; ++i){
                if(updated_level_set3[i] <= 0){
                    index_LG3 = i;
                }else{
                    W_vec_real[i] = W4_vec[i];
                }
            }
        }
        MATRIXreal.push_back(W_vec_real);
        MATRIX1.push_back(W1_vec_new);
        MATRIX2.push_back(W2_vec_new);
        if(num_states >= 3){MATRIX3.push_back(W3_vec_new);}
        if(num_states >= 4){MATRIX4.push_back(W4_vec_new);}
        W1_vec = W1_vec_new;
        W2_vec = W2_vec_new;
        if(num_states >= 3){W3_vec = W3_vec_new;}
        if(num_states >= 4){W4_vec = W4_vec_new;}
        level_set = updated_level_set;
        if(num_states >= 3){level_set2 = updated_level_set2;}
        if(num_states >= 4){level_set3 = updated_level_set3;}
        T += dt;        //time update
        time_vec.push_back(T);
        nT += 1;

    }

    cout << "(actual final time = " << T << "s)" << endl; 
    
    return MATRIXreal;
}


Vector writetofile(TestCase CASE, vector< vector<Prim> > MATRIX){
    Vector output(10); //returning: {u_min, u_max, rho_min, rho_max,
                                // p_min, p_max, e_min, e_max, c_min, c_max}
                                // for use in automatic plot scaling later
    Vector space = CASE.space;
    string name = CASE.datfile;
    int n = CASE.n;     //number of spatial nodes
    int nT = CASE.nT;   //number of time steps
    assert(MATRIX.size() == nT);
    assert(MATRIX[0].size() == n);
    ofstream ufile, rhofile, pfile, efile, gammafile, finalT;
    ufile.open("./data/u_" + name);
    rhofile.open("./data/rho_" + name);
    pfile.open("./data/p_" + name);
    efile.open("./data/e_" + name);
    gammafile.open("./data/gamma_" + name);
    finalT.open("./data/finalT" + name);
    double u_min=CASE.W1.u, u_max=u_min, rho_min=CASE.W1.rho, rho_max=rho_min; 
    double p_min=CASE.W1.p, p_max=p_min, e_min=p_min/(rho_min*(CASE.W1.gamma-1)), e_max=e_min;
    double gamma_min=CASE.W1.gamma, gamma_max=gamma_min;
    //just initialising these to possible values
    double e;
    for(int i = 0; i < n; ++i){
        //first column of each file is the space vector
        ufile << space[i] << ' ';
        rhofile << space[i] << ' ';
        pfile << space[i] << ' ';
        efile << space[i] << ' ';
        gammafile << space[i] << ' ';
        finalT << space[i] << ' ';
        for(int j = 0; j < nT; ++j){
            ufile << MATRIX[j][i].u << ' ';
            if(MATRIX[j][i].u < u_min){
                u_min = MATRIX[j][i].u;
            }
            if(MATRIX[j][i].u > u_max){
                u_max = MATRIX[j][i].u;
            }
            if(j== nT-1){
                finalT << MATRIX[j][i].u << ' ';
            }

            rhofile << MATRIX[j][i].rho << ' ';
            if(MATRIX[j][i].rho < rho_min){
                rho_min = MATRIX[j][i].rho;
            }
            if(MATRIX[j][i].rho > rho_max){
                rho_max = MATRIX[j][i].rho;
            }
            if(j== nT-1){
                finalT << MATRIX[j][i].rho << ' ';
            }

            pfile << MATRIX[j][i].p << ' ';
            if(MATRIX[j][i].p < p_min){
                p_min = MATRIX[j][i].p;
            }
            if(MATRIX[j][i].p > p_max){
                p_max = MATRIX[j][i].p;
            }
            if(j== nT-1){
                finalT << MATRIX[j][i].p << ' ';
            }

            e = MATRIX[j][i].p/(MATRIX[j][i].rho*(MATRIX[j][i].gamma-1));
            efile << e << ' ';
            if(e<e_min){
                e_min = e;
            }
            if(e>e_max){
                e_max = e;
            }
            if(j== nT-1){
                finalT << e << ' ';
            }
            gammafile << MATRIX[j][i].gamma << ' ';
            if(MATRIX[j][i].gamma < gamma_min){
                gamma_min = MATRIX[j][i].gamma;
            }
            if(MATRIX[j][i].gamma > gamma_max){
                gamma_max = MATRIX[j][i].gamma;
            }
            if(j== nT-1){
                finalT << MATRIX[j][i].gamma << ' ';
            }
            //essentially transposing matrix to file write s.t.
            //vector downwards in space, with each column stepping in time
        }
        ufile << '\n';
        rhofile << '\n';
        pfile << '\n';
        efile << '\n';
        gammafile << '\n';
        finalT << '\n';
    }
    ufile.close();
    rhofile.close();
    pfile.close();
    efile.close();
    gammafile.close();
    double tol = 1e-3;
    if(fabs(u_max - u_min) < tol){ u_min += -0.1; u_max += 0.1;}
    if(fabs(rho_max - rho_min) < tol){ rho_min += -0.1; rho_max += 0.1;}
    if(fabs(p_max - p_min) < tol){ p_min += -0.1; p_max += 0.1;}
    if(fabs(e_max - e_min) < tol){ e_min += -0.1; e_max += 0.1;}
    if(fabs(gamma_max - gamma_min) < tol){ gamma_min += -0.1; gamma_max += 0.1;}

    double arr[10] = {u_min, u_max, rho_min, rho_max, p_min, p_max, e_min, e_max, 
        gamma_min, gamma_max};
    output.assign(arr, &arr[10]);
    return output;
}

void writetoGNUplot(TestCase CASE, Vector domain){
    Vector space = CASE.space;
    string gpname = CASE.gpname;
    string datfile = CASE.datfile;
    string title = CASE.title;
    string giffile = CASE.giffile;
    int nT = CASE.nT;

    string u_datfile = "./data/u_" + datfile;
    string rho_datfile = "./data/rho_" + datfile;
    string p_datfile = "./data/p_" + datfile;
    string e_datfile = "./data/e_" + datfile;
    string gamma_datfile = "./data/gamma_" + datfile;

    ofstream file;
    file.open(gpname);
    double maxU = domain[1]; 
    double minU = domain[0];
    double maxX = *max_element(space.begin(), space.end());
    double minX = *min_element(space.begin(), space.end());
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal gif animate delay 10 \n";
    file << "n = " << nT << '\n';
    file << "set out \'" << giffile << "\' \n";
    file << "set key off \n";
    
    file << "set xrange[" << minX << ":" << maxX << "] \n";
    file << "do for [i=0:" << nT-1 << "] { \n    ";
    file << "j = i+2 \n";
    file << "set multiplot layout 2,2 \n";
    //rho plot:
    file << "set title \"" << "Density - rho" << "\" \n";
    file << "set yrange[" << domain[2] - 0.3*fabs(domain[3]-domain[2]) << \
        ":" << domain[3] + 0.3*fabs(domain[3]-domain[2]) << "] \n";
    file << "plot \"" << rho_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";
    //u plot:
    file << "set title \"" << "Velocity - u" << "\" \n";
    file << "set yrange[" << minU - 0.3*fabs(maxU-minU) << ":" << maxU + 0.3*fabs(maxU-minU) << "] \n";
    file << "plot \"" << u_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";
    /*
    //gamma plot:
    file << "set title \"" << "gamma" << "\" \n";
    file << "set yrange[" << domain[8] - 0.3*fabs(domain[9]-domain[8]) << \
        ":" << domain[9] + 0.3*fabs(domain[9]-domain[8]) << "] \n";
    file << "plot \"" << gamma_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";
    */
    //p plot:
    file << "set title \"" << "Pressure - p" << "\" \n";
    file << "set yrange[" << domain[4] - 0.3*fabs(domain[5]-domain[4]) << \
        ":" << domain[5] + 0.3*fabs(domain[5]-domain[4]) << "] \n";
    file << "plot \"" << p_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";
    //e plot:
    file << "set title \"" << "Internal energy - e" << "\" \n";
    file << "set yrange[" << domain[6] - 0.3*fabs(domain[7]-domain[6]) << \
        ":" << domain[7] + 0.3*fabs(domain[7]-domain[6]) << "] \n";
    file << "plot \"" << e_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";

    
    file << "unset multiplot \n } \n";
    file << "exit";
    file.close();
}

void GNUplot_finalT(TestCase CASE, string filetype){
    ofstream file;
    file.open(CASE.title + " finalT.gp");
    string datfile = "./data/finalT" + CASE.datfile;
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal " << filetype << " \n";
    file << "set xrange[0:1] \n";
    file << "set out './visualisation/"<< CASE.title << " finalT." << filetype <<"\' \n";
    file << "set key off \n";
    file << "set multiplot layout 2,2 \n";
    file << "set title \"Density - rho\" \n";
    file << "plot '"<< datfile << "\' using 1:3 with linespoints pt 6 ps 0.3 lc rgb \'blue\' \n";
    file << "set title \"Velocity - u\" \n";
    file << "plot '"<< datfile << "\' using 1:2 with linespoints pt 6 ps 0.3 lc rgb \'blue\' \n";
    file << "set title \"Pressure - p\" \n";
    file << "plot '"<< datfile << "\' using 1:4 with linespoints pt 6 ps 0.3 lc rgb \'blue\' \n"; 
    file << "set title \"Internal energy - e\" \n";
    file << "plot '"<< datfile << "\' using 1:5 with linespoints pt 6 ps 0.3 lc rgb \'blue\' \n";
    file << "exit";
}

void tellmethings(TestCase CASE, int nG, bool set_all_G, 
        bool re_init, bool MUSCL, bool isobar_fix){
    cout << "       ###, ,##, ,##,\n        #  # #  # #  #\n        ###  #  # #  #\n        #  # #  # #  #\n        ###' '##' '##'\n             .--,\n           /  (\n           /    \\ \n          /      \\ \n         /  0  0  \\ \n ((()   |    ()    |   ()))\n \\  ()  (  .____.  )  ()  /\n  |` \\_/ \\  `""`  / \\_/ `|\n  |       `.'--'.`       |\n   \\        `""`        /\n    \\                  /\n     `.              .'    ,\n      |`             |  _.'|\n      |              `-'  /\n      \\                 .'\n       `.____________.-' " << endl;
    cout << "|| " << CASE.title << " || " << endl;
    cout << "Basic Ghost fluid solver" << endl; 
    cout << "spatial discretisation = " << CASE.n << endl;
    cout << "simulated time = " << CASE.Tf << "s" << endl;
    cout << "number of ghost cells set to: ";
    if(set_all_G){ cout << "-all- \n";} else{cout << nG << endl;}
    cout << "Isobaric fix condition is: ";
    if(isobar_fix){ cout << "ON \n";}else{ cout<< "OFF \n";}
    cout << "level set reinitialisation set to: ";
    if(re_init){ cout << "ON \n";}else{ cout<< "OFF \n";}
    cout << "second order extension scheme is set to: ";
    if(MUSCL){ cout << "MUSCL with Superbee slope limiter \n";
    }else{ cout << "OFF - first order \n";} 
    cout << "number of time steps in simulated time = " << CASE.nT << endl;
    cout << "6 files generated in ./data/: \n  *" << CASE.datfile << "(x5 properties)\n";
    cout << "plus 1x all properties at final T." << endl; 
    cout << "for final time solution image in ./visualisation/: \n  *" << CASE.title << " finalT.svg" << endl;
    cout << "compile with: \"chmod u+x " << CASE.title << " finalT.gp\" then run." << endl;
    cout << "to generate a gif anaimation of full solution evolution in ./visualisation: \n";
    cout << "  *" << CASE.gpname << endl;
    cout << "compile with: \"chmod u+x " << CASE.gpname << \
        "\" then run. " << endl;
    cout << "done." << endl;
    
}


int main(){
    
    //~~ DEFINE TEST PARAMETERS HERE ~~//
    TestCase CASE("TestB");           //test case as defined in TestCases_GhostFluid.H
                                        //refer to .H file for string options
    //user specified settings:
    int nG = 3;                         //number of ghost cells
    bool set_all_G = 1;                 //override num ghost cells and set all cells: on/off 
    bool re_init = 1;                   //re-initialise level-set function at every iteration
    bool isobar_fix = 1;                //isobaric fix for ghost fluid method: on/off
    bool MUSCL = 1;                     //second order reconstruction with MUSCL: on/off
    string scheme = "SuperBee";         //slope limiter paired with MUSCL linear reconstruction
              

    
    Vector time_vec;
    vector< vector<Prim> > MATRIX = generate_solution(CASE, time_vec, nG, set_all_G, isobar_fix,
                                                        MUSCL, scheme, re_init);
    CASE.nT = MATRIX.size();
    
    //print simulation details to terminal:
    tellmethings(CASE, nG, set_all_G, re_init, MUSCL, isobar_fix);
    
    //write data and visualisation files:
    Vector domain = writetofile(CASE, MATRIX);
    writetoGNUplot(CASE, domain); 
    GNUplot_finalT(CASE, "svg");
    
    return 0;
}





