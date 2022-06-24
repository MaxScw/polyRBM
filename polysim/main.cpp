#include <iostream>
#include "r250.h"
#include "point.h"
#include "numerics.h"
#include <deque>
#include <map>
#include <sstream>
#include <fstream>

typedef IPoint V;

using namespace std;

template<> string V::coordSep(",");

# define LBIT   (8)
# define LSIZE  (1u<<LBIT)
# define NCELLS (LSIZE*LSIZE*LSIZE)
# define LMASK  (LSIZE-1)

static R250 r250;

class Latt
{
    char* cells;
    
public:
    
    Latt(){
        cells = new char[NCELLS];
        this->setZero();
    }
    
    ~Latt(){delete[] cells;}

    const char& operator[](const V& p) const
    {
        return cells[calcAddress(p)];
    }
    
    char& operator[](const V& p)
    {
        return cells[calcAddress(p)];
    }
    
    void setZero(){
        for ( int i = 0 ; i < NCELLS;++i) cells[i] = 0;
    }
    
    int calcAddress (const V& p) const
    {
        // cout << (p.X & LMASK) << " " << ((p.Y & LMASK)<<LBIT) << " " <<  ((p.Z & LMASK)<<(2*LBIT)) << endl;
        return (p.X & LMASK) | ((p.Y & LMASK)<<LBIT) | ((p.Z & LMASK)<<(2*LBIT));
    }
    
    V calcPosition (const int& address) const
    {
        return V( address & LMASK, (address >> LBIT ) & LMASK, (address >> 2*LBIT ) & LMASK );
    }
    
    void printOccupancy( ostream& strm = cout )
    {
        int count = 0 ; 
        for ( int i = 0 ; i < NCELLS; ++ i )
        {
           char val = cells[i];
           if ( val ) 
           {
               strm << "position " << calcPosition(i) << " --> " << int(val) << endl;
               count++;
           }
        }
        cout << "total = " << count << endl;
    }
};

struct Setup{
    int N,tt,dt,tmin; // chain length, max simulation time, minimal time for analysis
    bool is2D;
    bool isPhantom;
    bool printRing;
    bool notReturning;
    
    Setup(int argc,char *argv[])
    :N      ( argc > 1 ? atoi(argv[1]):32)
    ,is2D   ( argc > 2 ? atoi(argv[2]):false )// reduce to 2 dimensions?
    ,isPhantom( argc > 3 ? atoi(argv[3]):false )// ignore excluded volume?
    ,tt     ( argc > 4? atoi(argv[4]):10010000)
    ,dt     ( argc > 5 ? atoi(argv[5]):100)
    ,tmin   ( argc > 6 ? atoi(argv[6]):100000)
    ,printRing(argc>7? atoi(argv[7]):false )
    ,notReturning(argc>8? atoi(argv[8]):false )
    {
        if (argc<2)
        {
            cout << "usage : ./polysim  N  is2D  isPhantom  total_steps  interval_steps  relax_steps  print_ring_conformations not_returning_walk\n";
            exit(1);
        }
    }
    
    string genName (const string& pref = "", const string& suff = "") const
    {
        ostringstream strm;
        strm<<pref<<"N"<<N<<(is2D?"_2d":"_3d")<<(isPhantom?"_phantom":"_exvol")<<(printRing?"_ring":"_open")<<(notReturning?"_notreturn":"")<<suff;
        return strm.str();
    }
};

struct BondEncoder
{
    virtual string operator()(const V& b)=0;
};

struct BondToString : public BondEncoder 
{
    virtual string operator()(const V& b)
    {
        ostringstream o;
        o << b;
        return o.str();
    }
};

// struct BondToBits : public BondEncoder 
// {
//     virtual string operator()(const V& b)
//     {
//         
//         ostringstream o;
//         for ( int i = 0 ; i < 3-is2D; ++i )
//         {
//             o << b;
//         }
//         return o.str();
//     }
// };

struct Simulator
{
    const Setup     setup;
    Latt            latt;
    deque < V >     conf;
    size_t          time;
    
    typedef map <string,ofstream*> FMap;
    
    FMap files;
//     bool newFile;
    
    Simulator(const Setup& setup)
    :setup(setup)
    ,latt()
    ,time(0)
//     ,newFile(true)
    {}
    
    ~Simulator()
    {
        deleteFiles();
    }
    
    void deleteFiles()
    {
        for (FMap::iterator it = files.begin(); it != files.end(); ++it)
        {
            it->second->close();
            delete it->second;it->second=0;
        }
        files.clear();
    }
    
    void init()
    {
        conf.clear();
        latt.setZero();
        time = 0;
        
        this->deleteFiles();
//         newFile = true;
    }
    
    V getRandomBond()
    {
        V b;
        int dir=r250.Uniform()*(2+(!setup.is2D));
        int sgn=r250.Uniform()>0.5?1:-1;
        b[dir]=sgn;
        return b;
    }
    
    bool randomWalk()
    {           
        V pos(r250.Uniform()*LSIZE,r250.Uniform()*LSIZE,(setup.is2D)*r250.Uniform()*LSIZE);
        
        conf.push_back(pos);
        if ( !setup.isPhantom ) latt[pos]=1;
        
        int simuCount = 0;
        
        while ( conf.size() < setup.N && simuCount < 100 )
        {   
            if ( ! randomWalkStep() )
            {   
                this->simulate(100);
                ++simuCount;
            }
        }
        if (simuCount==100) return false;
        else return true;
    }

    bool wouldNotReturn(const V& pos) const
    {
        if ( conf.size() > 1 )
        {
            return !(conf[conf.size()-2] == pos);
        }
        else return true;
    }
    
    bool randomWalkStep()
    {
        int tryCount=0;
        while ( tryCount < 10 )
        {
            V pos  = conf.back();
            V bond = getRandomBond();
            pos   += bond;
            
            if ( (setup.isPhantom || latt[pos] == 0) && ( !setup.notReturning || wouldNotReturn(pos) ) )
            {
                conf.push_back(pos);
                if(!setup.isPhantom) latt[pos] = 1;
                return true;
            }
            ++tryCount;
        }
        return false;
    }
    
    void simulate(int n)
    {
        uint64_t accept=0;
        uint64_t nsweeps=n*uint64_t(conf.size());
        static bool frontMode = false;
        
        for ( int i = 0; i < nsweeps; ++i)
        {
            V p = frontMode ? conf.back() : conf.front() ;
            
            p  += this->getRandomBond();
            
            bool valid = true;
            
            if ( setup.notReturning && conf.size() > 1)
            {
                int returnIdx = (frontMode?conf.size()-2:1);
                if ( conf[returnIdx] == p ) valid = false;
            }
            
            if ( valid && !setup.isPhantom )
            {
                valid = latt[p] == 0;
            }

            if ( valid )
            {
                if ( frontMode ) 
                {
                    if(!setup.isPhantom)
                    {
                        latt[p] = 1;
                        latt[conf.front()] = 0;
                    }
                    conf.push_back(p);
                    conf.pop_front();
                }
                else 
                {
                    if(!setup.isPhantom)
                    {
                        latt[p] = 1;
                        latt[conf.back()] = 0;
                    }
                    conf.push_front(p);
                    conf.pop_back();
                }
                ++accept;
            }
            else
            {
                frontMode = !frontMode;
            }
        }
//         cout << " accept = " << accept <<" of " << nsweeps << endl;
        time+=n;
    }
    
    void printBonds(ostream& strm = cout, const string& sep=",", const string& lsep = ",", BondEncoder* enc = 0 ) const
    {
        V::setSeparator(sep);
//         strm << (lsep=="\n"?"#":"")<< time<<lsep;
        if ( lsep!="\n") strm << time << lsep;
        for ( int i = 0; i < conf.size()-1;++i)
        {
            if ( enc ) strm<<(*enc)(conf[i+1]-conf[i]);
            else strm<<(conf[i+1]-conf[i]);
            if (i<conf.size()-2)strm<<lsep;
        }
        strm<<endl;
    }
    
    void printConformation(ostream& strm = cout, const string& sep=",", const string& lsep = ",") const
    {
        V::setSeparator(sep);
        V first (conf.front());
        if ( lsep!="\n") strm << time << lsep;
        for ( int i = 0; i < conf.size();++i)
        {
            strm<<conf[i]-first;
            if (i<conf.size()-1)strm<<lsep;
        }
        strm<<endl;
    }
    
    ofstream* getFile(const string& name)
    {
        FMap::iterator it = files.find(name);
        if ( it == files.end() )
        {
            files[name] = new ofstream(name.c_str(),ios::out);
        }
        return files[name];
    }
    
    
    void writeFiles(const ios::openmode& mode=ios::app) 
    {
        ostringstream endstrm;
        endstrm<<"_";
        string lsep=",";
        string sep=",";
        
        if ( mode != ios::app )
        {
            lsep="\n";
            sep=" ";
            endstrm<<setw(7)<<setfill('0')<<time<<"_";
        }
        string endStr = endstrm.str();
        
//         ios::openmode mymode=mode;
// //         if ( newFile ) { mymode=ios::out; const_cast<Simulator*>(this)->newFile = false;}
   
        ofstream* bf(this->getFile(setup.genName("",endStr+"bonds.csv")));
        ofstream* cf(this->getFile(setup.genName("",endStr+"conf.csv")));
        
        this->printBonds(*bf,sep,lsep);
        this->printConformation(*cf,sep,lsep);
        
//         bf.close();
//         cf.close();
        
    }
    
    bool isRing() const
    {
        V diff = conf.back()-conf.front();
        int diffSqr(diff*diff);
        return (diffSqr==1);
    }
    
};

int main(int argc, char *argv[])
{
    
    Simulator sim(Setup(argc,argv));
    cout << "running setup " << sim.setup.genName()<<endl;
    
    for ( int i = 0; i < 10000; ++i ) r250.Refresh();
    sim.randomWalk();

    Correlator < Product<DPoint,double,double> , MultiplyIncrementor<> > reeCorr;
    M2_new<double>Rg;
    
    cout << "lbits = " << LBIT << " lmask = " << LMASK << " lsize = " << LSIZE << " NCELLS = " << NCELLS<< endl;
    
    int sampleCount=0;
    while ( sim.time <= sim.setup.tt )
    {
        sim.simulate(sim.setup.dt);
        if ( sim.time > sim.setup.tmin ) 
        {
            if ( sim.setup.printRing==false || sim.isRing() )
            {
                sim.writeFiles(ios::app);
                
                if ( sampleCount % 10000 == 0) 
                {
                    ofstream l(sim.setup.genName("","_lastconf.dat").c_str());
                    sim.printConformation(l," ","\n");
                }
                
                V Re=sim.conf.back()-sim.conf.front();
//             reeCorr(Re);
                DPoint com;
                com=CenterOfMass(sim.conf.begin(),sim.conf.end());
                Rg+=SqrRadiusOfGyration(sim.conf.begin(),sim.conf.end(),com);
            }
            if (sampleCount++ % 10000 == 0){
                cout << "time = " << sim.time << endl;
//                 sim.printConformation(cout," ","\n");
            }
        }
    }
//     string reeName = sim.setup.genName("","_reeCoor.dat");  
//     ofstream out1(reeName.c_str()); reeCorr.Print(out1);

    string rgName = sim.setup.genName("","_Rg2.dat"); ofstream out2(rgName.c_str()); 
    out2<<sim.conf.size()<<" "<<Rg.get1()<<" "<<Rg.get2()<<" "<<Rg.get0()<<endl;

    return 0;
}

