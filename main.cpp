#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

//myElements[i].wezly niepotrzebne??!?!?!?!?!

const int ne = 12;//elements
const int nh = 21;//wezly, 0-20

double Hglob[nh][nh] = {};	//globalna macierz wspolczynnikow ukladow rownan
double Pglob[nh] = {};		//prawa czesc ukladu rownan
double temperatures[nh] = {};			//wektor niewiadomych (temperatur)

double eps = 1e-12;

/////////soe/////////
double systemOfEquations[nh][nh];

void setSystemOfEquations() {

	for(int i=0; i<nh+1; i++)
		for(int j=0; j<nh; j++)
		{
			if(i != nh)
				systemOfEquations[j][i] = Hglob[j][i];
			else
				systemOfEquations[j][i] = Pglob[j];
		}
}

// nie mam pewnosci ze dobrze
bool solveSOE() {
	int max;
	double tmp;
	int N = nh;
	double a[nh][nh];

	for(int i=0; i<nh; i++)
		for(int j=0; j<nh; j++)
			a[i][j] = Hglob[i][j];

	double m, s;

	for(int i=0; i<nh - 1; i++) // eliminacja
	{
		for(int j=i+1; j<N; j++)
		{
			if(fabs(a[j][i]) < eps) return false;
				m = -a[j][i] / a[i][i];	
			for(int k=i + 1; k<= N; k++)
				a[j][k] += m * a[i][k];
		}
	}

	for(int i = nh - 1; i >= 0; i--)
	{ 	
		s = a[i][nh];

		for(int j=N-1; j>=0; j--)
		{
			for(int k = nh - 1; k >= i + 1; k--) s -= a[i][k] * temperatures[k]; 
				if(fabs(a[i][i]) < eps)
					return false; 
				temperatures[j] = s/a[i][i];
		}
	}
	return true;
}
///////po soe////////

//grid &grid1,double **Hg,double *Pg,node *tab1
void gauss()
{

    for (int k = 0; k<ne; k++){   //gaus -
            double dziel = Hglob[k][k];
            for (int i = (k+1); i<nh; i++)
            {
            	double wsp = Hglob[i][k] / Hglob[k][k];

            	for (int j = k; j<nh; j++) // idzie w prawo
            		Hglob[i][j] -= Hglob[k][j]*wsp;
            	Pglob[i] -= Pglob[k] * wsp;
            }
    }

    for (int i = ne; i >= 0; i--){   //obliczanie temp od dołu macierzy
            double suma = 0;
            for (int j = (i+1); j<nh; j++)//końcowy gaus
                    suma = suma + Hglob[i][j] * temperatures[i];
            temperatures[i] = (-Pglob[i] - suma) / Hglob[i][i];
    }
}

// argumenty: Hglob, nh
void Gauss()//double **macierz, int n
{
    double s, m;

    for(int i=0; i < ne-1; i++)
            for(int j = i+1; j < ne; j++)
            {
                    m = Hglob[j][i] / Hglob[i][i];

                    for(int k=i+1; k <= ne; k++)
                            Hglob[j][k] += m * Hglob[i][k];
            }

    double *tab=new double[ne];
   
    for(int i=0; i<ne; i++)
            tab[i] = 0;

    for(int i = ne-1; i>=0; i--)
    {
            s = Hglob[i][ne];

            for(int j=ne-1; j>=i+1; j--)
                    s -= Hglob[i][j]*tab[j];

            tab[i] = s / Hglob[i][i];
    }
    
    //for(int i=0; i<ne; i++)
               //for(int j=0; j<ne; j++)
                       //cout << Hglob[i][j] << "\t";
    for(int i=0; i<ne; i++)
    {
    	cout << tab[i] << " ";
    }

}

struct node {
	//wsp. lokalne:
	int x;
	int y;

	//globalny nr wezla
	int index;
};

class Element
{
      public:
             node nodes[4];		// names of nodes
             int typeOfBC;		// 0 or 1 or 2
             int whereBC;		// 0 or 1 or 2 or 3 or 4

             double H[4][4];		// macierz lokalna
             double P[4];			//lokalna prawej czesci ukladu rownan
} myElements[12] = {}; // tablica globalna, z³o |..|

// shows all elements and their data
void showMyElements()
{
     for(int i=0; i<ne; i++)
     {
              cout << "element nr " << i << ": " << endl;

              cout << "wspolrzedne globalne dla wierzcholkow:" << endl;
              for(int j=0; j<4; j++)
              	cout << j << ": (" << myElements[i].nodes[j].x << ", " << myElements[i].nodes[j].y << ")" << endl;

              cout << "wsp. glob. dla wierzcholkow [2]:" << endl;
              for(int j=0; j<4; j++)
              	cout << myElements[i].nodes[j].index << " ";

              cout << endl;

              cout << endl << "Typ BC: " << myElements[i].typeOfBC << " na sciance: " << myElements[i].whereBC;

              cout<<endl<<"macierz lokalna " << i << ":" << endl;
              for(int j=0; j<4; j++){
              	for(int k=0; k<4; k++)
              		cout << " " << myElements[i].H[j][k];
              	cout << endl;
              }
              cout << "P lokalne:" << endl;
              for(int j=0; j<4; j++)
              	cout << myElements[i].P[j] << " ";
              cout << endl << endl;
    }
}

void showTemperatures()
{
	cout << endl << endl << "TEMPERATURY:" << endl;
	for(int i=0; i<nh; i++)
		cout << temperatures[i] << " ";
	cout << endl;
} 

//double dndxl = pochFunKsztaltu(l, ksi[j][0], ksi[j][1], 1) * Jodwr[0][0] + pochFunKsztaltu(l, ksi[j][0], ksi[j][1], 2) * Jodwr[0][1];
double pochFunKsztaltu(int n, double ksi, double ni, int type)
{
       double result;
       if(type == 1) {
               switch(n) {
                         case 0:
                              result = -0.25 * (1.0 - ksi);
                              break;
                         case 1:
                              result = 0.25 * (1.0 - ksi);
                              break;
                         case 2:
                              result = 0.25 * (1.0 + ksi);
                              break;
                         case 3:
                              result = -0.25 * (1.0 + ksi);
                              break;
                         }
               }
       else {
            switch(n) {          
                case 0:
                     result = -0.25 * (1.0 - ni);
                     break;
                case 1:
                     result = -0.25 * (1.0 + ni);
                     break;
                case 2:
                     result = 0.25 * (1.0 + ni);
                     break;
                case 3:
                     result = 0.25 * (1.0 - ni);
                     break;
                 }
            }
            
       return result;
}

// loads data from dane.txt
void loadFromFile()
{
     int tmp;
     
     ifstream plik;
     plik.open("dane.txt"); // zakladam ze nie ma erroru
     
     for(int i=0; i<ne; i++)
     {
              plik >> myElements[i].nodes[0].x;
              plik >> myElements[i].nodes[0].y;
              plik >> myElements[i].nodes[1].x;
              plik >> myElements[i].nodes[1].y;
              plik >> myElements[i].nodes[2].x;
              plik >> myElements[i].nodes[2].y;
              plik >> myElements[i].nodes[3].x;
              plik >> myElements[i].nodes[3].y;

              for(int j=0; j<4; j++)
              	plik >> myElements[i].nodes[j].index;

              plik >> myElements[i].typeOfBC;
              plik >> myElements[i].whereBC;
    }
    plik.close();
}

void showHglob() {

	cout << "MACIERZ GLOBALNA:" << endl;
	for(int i=0; i<nh; i++)
	{
		for(int j=0; j<nh; j++)
			cout << setprecision(4) << Hglob[i][j] << "\t";
			
		cout << endl;
	}
	cout << endl;

	cout << "P globalne: " << endl;
	for(int i=0; i<nh; i++)
		cout << Pglob[i] << " ";
}

int main()
{
    /* -- zmienne poczatkowe -- */
    int t0 = 100;        // temperatura poczatkowa          [st. C]
    int tOt = 1200;      // temperatura otoczenia           [st. C]
    unsigned int K = 25; // wsp. przewodzenia ciepla        [mK]
    double alfa = 300; // wsp. wymiany ciepla             [W / m^2]
    unsigned q = 1000;   // strumien ciepla                 [W]
    
    // Macierz Jacobiego – macierz zbudowana z pochodnych czastkowych (pierwszego rzedu) funkcji, której skladowymi sa funkcje rzeczywiste
    double J[2][2] = {}; // zerowanie
    // Macierz odwrotna
    double Jodwr[2][2];

    double L = 1;
    double tmpuuu[2][2];//wzor 46
	tmpuuu[0][0] = 100;//2/6 * alfa * L;//2/6 * 300 == 100
  	tmpuuu[0][1] = 50;//alfa * L/6;//1/6 * 300 == 50
  	tmpuuu[1][0] = 50;//alfa * L / 6;//50
  	tmpuuu[1][1] = 100;//2/6 * alfa * L;//100
    
    double ksi[4][2] = {// wsp. pktow calkowania Gaussa dla calek powierzchniowych, ksi ale TAKZE ni
           { -0.5773502692,  -0.5773502692 },
           {  0.5773502692,  -0.5773502692 },
           {  0.5773502692,   0.5773502692 },
           { -0.5773502692,   0.5773502692 }
    };
    loadFromFile();// wczytanie danych z pliku

    for(int i=0; i<ne; i++)// main petla po elementach!!!!!!!!!!!!!!!!!!!11111111111111111111111111111111
    {
		for(int j=0; j<4; j++)// zerowanie macierzy dla kazdego elementu
		{
			for(int k=0; k<4; k++)
				myElements[i].H[j][k] = 0;//zerowanie H lok

			myElements[i].P[j] = 0;//zerowanie P lok
		}

		for(int j=0; j<4; j++)//petla po pktach calkowania
		{
			for(int k=0; k<4; k++) //petla po wezlach, wzor 39: jakobian to macierz 2*2 sum od 0 do 3
			{
				//poch. f. ksztaltu - wzory 27 - 38
				// reszta: wzor 39!
				J[0][0] += myElements[i].nodes[k].x * pochFunKsztaltu(k, ksi[j][0], ksi[j][1], 1);
				J[0][1] += myElements[i].nodes[k].x * pochFunKsztaltu(k, ksi[j][0], ksi[j][1], 2);
				J[1][0] += myElements[i].nodes[k].y * pochFunKsztaltu(k, ksi[j][0], ksi[j][1], 1);
				J[1][1] += myElements[i].nodes[k].y * pochFunKsztaltu(k, ksi[j][0], ksi[j][1], 2);
			}//koniec petli po wezlach
			
			// wyznacznik potrzebny do wzoru 41
			double detJ = J[0][0] * J[1][1] - J[1][0] * J[0][1];
			//cout << detJ << endl;

			// obliczenie odwrotnego jakobianu, wzor 41
			Jodwr[0][0] = (1/detJ) * J[1][1];
			Jodwr[0][1] = (1/detJ) * (-J[0][1]);
			Jodwr[1][0] = (1/detJ) * (-J[1][0]);
			Jodwr[1][1] = (1/detJ) * J[0][0];
		//0.5, 0, 0, 0.5

		double Ndx[4], Ndy[4]; // wzor 42
			
		for(int k=0; k < 4; k++) // wzor 42, pochodne potrzebne do wzoru 43
		{
			Ndx[k] = pochFunKsztaltu(k, ksi[j][0], ksi[j][1], 1) * Jodwr[0][0] + pochFunKsztaltu(k, ksi[j][0], ksi[j][1], 2) * Jodwr[0][1];
			Ndy[k] = pochFunKsztaltu(k, ksi[j][0], ksi[j][1], 1) * Jodwr[1][0] + pochFunKsztaltu(k, ksi[j][0], ksi[j][1], 2) * Jodwr[1][1];
		}

		// double pochFunKsztaltu(unsigned int n, double ksi, double ni, int type)
		// typ 1: po ksi. typ 2: po ni
		for(int k=0; k<4; k++)
		   for(int n=0; n<4; n++)//poprawic wyswietlanie
		           myElements[i].H[k][n] += K * (Ndx[k] * Ndx[n] + Ndy[k] * Ndy[n]) * detJ; // wzor 43!
		       
		}//koniec petli po pktach calkowania
	}//koniec petli po elementach

      	// teraz 'tylko' wstawic macierze lokalne w odp. pozycje mac. globalnych i rozwiazac ukł. =ań
        
        //Ostatni etap rozw. zadania:
        //wstawienie macierzy lokalnych w odpowiednie pozycje w macierzy globalnej i rozwiązanie ukladu równań.
        //agregacja

        ///////////////war. brzegowe wstawiane do P lok//////////////////
		for(int i=0; i<ne; i++)//po elementach
		{
			//wzor 47!!!
			if(myElements[i].typeOfBC == 1)//strumien ciepla, dla gornej krawedzi
				if(myElements[i].whereBC == 3)//na trzeciej==gornej sciance elementu (czerwony)
				{
					cout << "wywoluje BC 1 dla i == " << i << endl;

					//dla wezlow prawy-gorny i lewy-gorny
					myElements[i].P[2] += 0.5 * q * L;//wzor 48!!!
					//cout << "dodane. obecnie wartosc to: " << myElements[i].P[3] << endl;
					myElements[i].P[3] += 0.5 * q * L;//wzor 48!!!
					//cout << "dodane. obecnie wartosc to: " << myElements[i].P[4] << endl;

					int nop1 = myElements[i].nodes[2].index;
					int nop2 = myElements[i].nodes[3].index;

					Hglob[nop1][nop1] += tmpuuu[0][0];
					Hglob[nop1][nop2] += tmpuuu[0][1];
					Hglob[nop2][nop1] += tmpuuu[1][0];
					Hglob[nop2][nop2] += tmpuuu[1][1];

				}
				
			if(myElements[i].typeOfBC == 2)//wymiana ciepla z otoczeniem
				if(myElements[i].whereBC == 2)// na drugiej==prawej sciance elementu (zielony)
				{
					cout << "BC dla i == " << i << endl;

					//dla wezlow prawy-dolny i prawy-gorny
					myElements[i].P[1] = 0.5 * alfa * tOt * L;// 0.5 * alfa * 1200 * 2 // 0.5 * 300 * 1200 * 1 == 180k
					myElements[i].P[2] = 0.5 * alfa * tOt * L;// 0.5 * alfa * 1200 * 2) // 0.5 * 300 * 1200 * 1 == 180k

					int nop1 = myElements[i].nodes[1].index;
					int nop2 = myElements[i].nodes[2].index;

					Hglob[nop1][nop1] += tmpuuu[0][0];
					Hglob[nop1][nop2] += tmpuuu[0][1];
					Hglob[nop2][nop1] += tmpuuu[1][0];
					Hglob[nop2][nop2] += tmpuuu[1][1];

				}
		}

		/////////////lokalne do globalnych///////////////////////////
		for(int i=0; i<ne; i++)//po wszystkich elementach
		{
			for(int j=0; j<4; j++)//po wezlach x
			{
				int tmp = myElements[i].nodes[j].index; //P globne
				Pglob[tmp] += myElements[i].P[j];		// P globalne

				for(int k=0; k<4; k++)//po wezlach y
				{
					int nop1 = myElements[i].nodes[j].index;
					int nop2 = myElements[i].nodes[k].index;

					Hglob[nop1][nop2] += myElements[i].H[j][k];

					//tak bylo ale zamiast wartosci trzeba brac numeracje
					//Hglob[i+j][i+k] += tmp;//myElements[i].H[j][k];
				}
			}
		}

        showMyElements();
        showHglob();

        for(int i=0; i<nh; i++)
        	temperatures[i] = 100;

		showTemperatures();

        setSystemOfEquations();
        solveSOE();
        //cout << endl << "po" << endl;
        //Gauss();
        //gauss();
        showTemperatures();

        // WYNIK: wektor temperatur == temperatury w kazdym wezle
        // wyniki: 1402 - 1203
        // http://home.agh.edu.pl/~byrska/IS/wynikiStacj.png
        //cin >> a;
}
//16, -4,16  0, 0, 0, -4,16
//-4,16

/*
taki powinny być wyniki P glob
0	0	0	0	500
0	0	0	0	1000
0	0	0	0	500
0	0	0	180000	360000
180000
*/

//100 100 100 100 -33.46 99.91 99.56 96.87 76.87 -18.42 96.61 93.78 83.38 39.94 -1.145 91.18 84.89 78.83 -1616 -1438 -1180 
//100 100 100 100 -33.46 99.91 99.56 96.87 76.87 -18.42 96.61 93.78 83.38 39.94 -1.145 91.18 84.89 78.83 -1616 -1438 -1180