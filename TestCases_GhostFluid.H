//Module: Toro Test input cases


class TestCase{
    public:
        int n;          //spatial discretisation- number of points
                        //1 boundary ghost cells are added to each side
        double L;       //1D spatial Length
        double Tf;      //final time
        double T0;      //initial time
        double CFL;     //Courant number
        double dx;      //spatial step
        double x0;      //property discontinuity position
        double x1;      //second discontinuity for 3 state cases
        double x2;      //third discontinuity for 4 states
        Prim W1;        //left side initialisation
        Prim W2;        //right side initialisation
        Prim W3;        //for cases with a third state
        Prim W4;        //for cases with a 4th state
        int num_states; //number of initial states ie. 2 or 3
        Vector space;   
        Vector level_set;
        Vector level_set2;  //second material interface needs to be tracked for 3 state cases
        Vector level_set3;  //for interface between states 3-4
        string gpname, datfile, title, giffile;
        int nT;         //number of time steps,
                        //not assigned until after solution and number of time steps are known
        TestCase(string CASE):
            W1(0,0,0,0), 
            W2(0,0,0,0),
            W3(0,0,0,0),
            W4(0,0,0,0) {
                if(CASE == "GhostTest1"){
                    title = "Test Ghost1";
                    datfile = "euler_eqs_ghosttest1.dat";
                    gpname = "gifplot_ghosttest1.gp";
                    giffile = "./visualisation/ghosttest1.gif";
                    num_states = 2;
                    W1.rho = 1.0; W1.u = 0.5; W1.p = 1.0; W1.gamma = 1.4;
                    W2.rho = 0.138; W2.u = 0.5; W2.p = 1.0; W2.gamma = 1.4;
                    n = 102;    //100 cells + 2 boundary ghost cells on each side
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 1.0;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2);
                    x0 = 0.25;
                }
                if(CASE == "GhostTest2"){
                    title = "Test Ghost2";
                    datfile = "euler_eqs_ghosttest2.dat";
                    gpname = "gifplot_ghosttest2.gp";
                    giffile = "./visualisation/ghosttest2.gif";
                    num_states = 2;
                    W1.rho = 1.0; W1.u = 0.5; W1.p = 1.0; W1.gamma = 1.4;
                    W2.rho = 0.138; W2.u = 0.5; W2.p = 1.0; W2.gamma = 1.67;
                    n = 102;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 1.0;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2);
                    x0 = 0.25;
                }
                if(CASE == "TestFig4"){
                    title = "Test Figure 4";
                    datfile = "euler_eqs_testfig4.dat";
                    gpname = "gifplot_testfig4.gp";
                    giffile = "./visualisation/testfig4.gif";
                    num_states = 2;
                    W1.rho = 2.0; W1.u = 0.0; W1.p = 9.8e5; W1.gamma = 1.4;
                    W2.rho = 1.0; W2.u = 0.0; W2.p = 2.45e5; W2.gamma = 1.4;
                    n = 102;
                    L = 4.0;    //space domain 0 < x < 1
                    Tf = 0.0022;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2);
                    x0 = 2.0;
                }
                if(CASE == "Toro411"){
                    //Toro tests in table 4.1 section 4.3.3
                    title = "Toro411";
                    datfile = "euler_eqs_Toro_411.dat";
                    gpname = "gifplot_Toro411.gp";
                    giffile = "./visualisation/Toro411.gif";
                    num_states = 2;
                    W1.rho = 1.0; W1.u = 0.0; W1.p = 1.0; W1.gamma = 1.4;
                    W2.rho = 0.125; W2.u = 0.0; W2.p = 0.1; W2.gamma = 1.4;
                    n = 102;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.25;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2);
                    x0 = 0.5;
                }
                if(CASE == "Toro412"){
                    //Toro tests in table 4.1 section 4.3.3
                    title = "Toro412";
                    datfile = "euler_eqs_Toro_412.dat";
                    gpname = "gifplot_Toro412.gp";
                    giffile = "./visualisation/Toro412.gif";
                    num_states = 2;
                    W1.rho = 1.0; W1.u = -2.0; W1.p = 0.4; W1.gamma = 1.4;
                    W2.rho = 1.0; W2.u = 2.0; W2.p = 0.4; W2.gamma = 1.4;
                    n = 102;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.15;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2);
                    x0 = 0.5;
                }
                if(CASE == "Toro413"){
                    //Toro tests in table 4.1 section 4.3.3
                    title = "Toro413";
                    datfile = "euler_eqs_Toro_413.dat";
                    gpname = "gifplot_Toro413.gp";
                    giffile = "./visualisation/Toro413.gif";
                    num_states = 2;
                    W1.rho = 1.0; W1.u = 0.0; W1.p = 1000.0; W1.gamma = 1.4;
                    W2.rho = 1.0; W2.u = 0.0; W2.p = 0.01; W2.gamma = 1.4;
                    n = 102;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.012;
                    T0 = 0;
                    CFL = 0.3;
                    dx = L/(n-1-2);
                    x0 = 0.5;
                }
                if(CASE == "TestA"){
                    title = "Test A";
                    datfile = "euler_eqs_testA.dat";
                    gpname = "gifplot_testA.gp";
                    giffile = "./visualisation/testA.gif";
                    num_states = 2;
                    W1.rho = 1.0; W1.u = 0.0; W1.p = 1.0e5; W1.gamma = 1.4;
                    W2.rho = 0.125; W2.u = 0.0; W2.p = 1.0e4; W2.gamma = 1.2;
                    n = 402;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 7.0e-4;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2);
                    x0 = 0.5;
                }
                if(CASE == "TestB"){
                    title = "TestB";
                    datfile = "euler_eqs_testB.dat";
                    gpname = "gifplot_testB.gp";
                    giffile = "./visualisation/testB.gif";
                    num_states = 3;
                    W1.rho = 1.3333; W1.u = 0.3535*pow(1e5,0.5); W1.p = 1.5e5; W1.gamma = 1.4;
                    W2.rho = 1.0; W2.u = 0.0; W2.p = 1.0e5; W2.gamma = 1.4;
                    W3.rho = 0.1379; W3.u = 0.0; W3.p = 1.0e5; W3.gamma = 1.67;
                    n = 202;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.0012;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2);
                    x0 = 0.05;
                    x1 = 0.50;
                }
                if(CASE == "TestB2"){
                    title = "TestB2";
                    datfile = "euler_eqs_testB2.dat";
                    gpname = "gifplot_testB2.gp";
                    giffile = "./visualisation/testB2.gif";
                    num_states = 4;
                    W1.rho = 1.3333; W1.u = 0.3535*pow(1e5,0.5); W1.p = 1.5e5; W1.gamma = 1.4;
                    W2.rho = 1.0; W2.u = 0.0; W2.p = 1.0e5; W2.gamma = 1.4;
                    W3.rho = 0.1379; W3.u = 0.0; W3.p = 1.0e5; W3.gamma = 1.67;
                    W4.rho = 1.0; W4.u = 0.0; W4.p = 1.0e5; W4.gamma = 1.4;
                    n = 202;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.0014;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2);
                    x0 = 0.05;
                    x1 = 0.40;
                    x2 = 0.60;
                }
                if(CASE == "TestB2R"){
                    title = "TestB2R";
                    datfile = "euler_eqs_testB2R.dat";
                    gpname = "gifplot_testB2R.gp";
                    giffile = "./visualisation/testB2R.gif";
                    num_states = 4;
                    W1.rho = 5.71429; W1.u = 3086.87; W1.p = 1.165e7; W1.gamma = 1.4;
                    W2.rho = 1.0; W2.u = 0.0; W2.p = 1.0e5; W2.gamma = 1.4;
                    W3.rho = 0.1379; W3.u = 0.0; W3.p = 1.0e5; W3.gamma = 1.67;
                    W4.rho = 1.0; W4.u = 0.0; W4.p = 1.0e5; W4.gamma = 1.4;
                    n = 402;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.0002;
                    T0 = 0;
                    CFL = 0.8;
                    dx = L/(n-1-2);
                    x0 = 0.05;
                    x1 = 0.40;
                    x2 = 0.60;
                }


                space.resize(n);
                level_set.resize(n);
                level_set2.resize(n);
                level_set3.resize(n);
                
                for(int i = 0; i<n; ++i){
                    space[i] = i*dx - dx;
                    level_set[i] = space[i] - x0;
                    level_set2[i] = space[i] - x1;
                    level_set3[i] = space[i] - x2;
                }

            }
};

