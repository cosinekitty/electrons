/*
    fastsim.cpp  -  Don Cross  -  1 July 2016.
    
    A simulation of mutually repulsive particles trapped on the surface of a sphere.
    This is also known as the Thompson Problem.
    See:  https://en.wikipedia.org/wiki/Thomson_problem
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
#include <cstdint>
#include <cmath>
#include <cstring>

namespace Electrons
{
    // Vector ---------------------------------------------------------------

    struct Vector
    {
        double  x;
        double  y;
        double  z;
        
        Vector(): x(0), y(0), z(0) {}
        Vector(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
        
        void Reset()
        {
            x = 0;
            y = 0;
            z = 0;
        }
        
        double MagSquared() const 
        {
            return (x*x) + (y*y) + (z*z);
        }
        
        double Mag() const
        {
            return sqrt(MagSquared());
        }
        
        Vector UnitVector() const
        {
            double r = Mag();
            if (r < 1.0e-9)
            {
                throw "Cannot find unit vector for near-zero magnitude vector.";
            }
            return Vector(x/r, y/r, z/r);
        }
        
        Vector& operator += (const Vector& other)
        {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }
        
        Vector& operator -= (const Vector& other)
        {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }
        
        static double Dot(const Vector& a, const Vector& b)
        {
            return (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
        }
    };
    
    typedef std::vector<Vector> VectorList;
    
    inline Vector operator + (const Vector& a, const Vector& b)
    {
        return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
    }
    
    inline Vector operator - (const Vector& a, const Vector& b)
    {
        return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
    }
    
    inline Vector operator * (double s, const Vector& v)
    {
        return Vector(s*v.x, s*v.y, s*v.z);
    }
    
    // Particle ---------------------------------------------------------------
    
    struct Particle
    {
        Vector position;
        Vector force;
        
        Particle(): position(), force() {}
        Particle(const Vector& _position): position(_position), force() {}
        Particle(const Vector& _position, const Vector& _force): position(_position), force(_force) {}
    };
    
    typedef std::vector<Particle> ParticleList;

    // JSON output ------------------------------------------------------------
    
    void Print(std::ostream& output, double t)
    {
        using namespace std;
        output << setprecision(14) << fixed;
        if (t >= 0.0) output << " ";
        output << t;
    }
    
    void JsonIndent(std::ostream& output, int indent)
    {
        int spaces = 4*indent;
        for (int i=0; i < spaces; ++i)
        {
            output << " ";
        }
    }
    
    void JsonPrint(std::ostream& output, const Vector& v)
    {
        output << "{\"x\":"; 
        Print(output, v.x); 
        output << ", \"y\":"; 
        Print(output, v.y);
        output << ", \"z\":";
        Print(output, v.z);
        output << "}";
    }
    
    void JsonPrint(std::ostream& output, const Particle& p)
    {
        output << "{\"position\":";
        JsonPrint(output, p.position);
        output << ", \"force\":";
        JsonPrint(output, p.force);
        output << "}";
    }
    
    void JsonPrint(std::ostream& output, const VectorList& list, int indent)
    {
        using namespace std;
        JsonIndent(output, indent);
        output << "[\n";
        bool first = true;
        for (const Vector& v : list) 
        {
            JsonIndent(output, indent);
            output << (first ? "    " : ",   ");
            JsonPrint(output, v);
            output << "\n";
            first = false;
        }
        JsonIndent(output, indent);
        output << "]\n";
    }
    
    // Random generators ------------------------------------------------------------
    
    double Random(std::ifstream& infile)
    {
        using namespace std;
        
        uint64_t data;
        infile.read((char *)&data, sizeof(data));
        if (!infile) throw "Read failure generating random number.";
        
        // Convert the 64-bit integer to double-precision floating point.
        // Divide by maximum possible value to get a number in the closed range [0, +1].
        double x = static_cast<double>(data) / static_cast<double>(UINT64_MAX);
        
        // Convert to the closed range [-1, +1].        
        return 1.0 - (2.0 * x);
    }

    Vector RandomSpherePoint(std::ifstream& infile)
    {    
        // Algorithm for picking a random point on a sphere.
        // Avoids any clustering of points.
        // See equations (9), (10), (11) in:
        // http://mathworld.wolfram.com/SpherePointPicking.html
        while (true)
        {
            double a = Random(infile);
            double b = Random(infile);
            double mag = (a*a) + (b*b);
            if (mag < 1.0)
            {
                double root = 2.0 * sqrt(1.0 - mag);
                return Vector(a*root, b*root, 1.0 - (2.0*mag));
            }
        }
    }

    // Simulation ---------------------------------------------------------------

    class Simulation
    {
    private:
        ParticleList particles;
        int frame;        
    
    public:
        Simulation(int numPoints)       // create an initial random state
            : frame(0)
        {
            using namespace std;
            
            if (numPoints < 0) throw "Number of points not allowed to be negative.";        
            ifstream infile("/dev/urandom", ios::in | ios::binary);
            if (!infile) throw "Could not open /dev/urandom to obtain random numbers.";
            particles.reserve(static_cast<ParticleList::size_type>(numPoints));
            for (int i=0; i < numPoints; ++i)
            {
                particles.push_back(Particle(RandomSpherePoint(infile)));
            }
        }
        
        int ParticleCount() const
        {
            return static_cast<int>(particles.size());
        }
        
        int FrameNumber() const
        {
            return frame;
        }
        
        void JsonPrint(std::ostream& output, int indent) const
        {
            using namespace std;
            
            JsonIndent(output, indent);
            output << "{\"frame\": " << frame << ", \"particles\": [\n";
            bool first = true;
            for (const Particle& p : particles) 
            {
                JsonIndent(output, indent);
                output << (first ? "    " : ",   ");
                Electrons::JsonPrint(output, p);
                output << "\n";
                first = false;
            }
            JsonIndent(output, indent);
            output << "]}\n";
        }
        
        bool Converge(double dt, int maxFrames)
        {
            const double forceTolerance = 1.0e-10;
            double force = 1;
            while (force > forceTolerance)
            {
                ++frame;
                CalcTangentialForces(particles);
                force = UpdatePositions(particles, particles, dt);
                
                //std::cout << "i=" << i << ", F=" << force << ", dt=" << dt << std::endl;
                if (frame > maxFrames)
                {
                    // Simulation did not converge within specified number of frames.
                    return false;
                }
            }
            return true;
        }
        
        bool AutoConverge()
        {
            // Theory #1: the simulation converges if total potential energy always decreases.
            // Theory #2: the simulation converges if max (particle tangential force) always decreases.
            
            
            
            return false;
        }
        
    private:
        static void CalcTangentialForces(ParticleList& particles)
        {
            // Reset each particle's force to a zero vector.
            for (Particle& p : particles)
            {
                p.force.Reset();
            }
            
            // Compute force between each unique pair of particles.
            const ParticleList::size_type numParticles = particles.size();
            for (ParticleList::size_type i=0; i < numParticles-1; ++i)
            {
                for (ParticleList::size_type j=i+1; j < numParticles; ++j)
                {
                    AddForces(particles[i], particles[j]);
                }
            }
            
            // The particles can only move along the surface of the sphere,
            // so eliminate the radial component of all forces, 
            // leaving tangential forces.
            for (Particle& p : particles)
            {
                // Calculate radial component using dot product and subtract
                // to get tangential component.
                p.force -= Vector::Dot(p.force, p.position) * p.position;
            }
        }
        
        static double UpdatePositions(ParticleList& inlist, ParticleList& outlist, double dt)
        {
            double maxForceMag = 0.0;
            const ParticleList::size_type numParticles = inlist.size();
            if (numParticles != outlist.size())
            {
                throw "Inconsistent particle list sizes.";
            }
            
            for (ParticleList::size_type i=0; i < numParticles; ++i)
            {
                Particle& p = inlist[i];
                Particle& q = outlist[i];
                
                double forceMag = p.force.Mag();
                if (forceMag > maxForceMag)
                {
                    maxForceMag = forceMag;
                }
                q.position = ((dt * p.force) + p.position).UnitVector();
            }
            return maxForceMag;
        }
        
        static void AddForces(Particle& a, Particle& b)
        {
            // Force of electrically charged particles:
            // F = k*q1*q2/r^2.
            // We simplify the problem as:  F = 1/r^2.
            // Force is along the direction of the line passing through both.
            Vector dp = a.position - b.position;
            double forcemag = 1.0 / dp.MagSquared();
            Vector force = forcemag * dp.UnitVector();
            a.force += force;
            b.force -= force;
        }
    };
}

//======================================================================================

class StandardDeviationCalculator
{
    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm    
private:
    int n;
    double mean;
    double sum;
    double max;
    double min;

public:
    StandardDeviationCalculator(): n(0), mean(0.0), sum(0.0), max(0.0), min(0.0) {}
    
    int NumData() const
    {
        return n;
    }
        
    void Append(double x)
    {
        ++n;
        if (n == 1)
        {
            max = min = x;
        }
        else
        {
            if (x < min) min = x;
            if (x > max) max = x;
        }
        double delta = x - mean;
        mean += delta/n;
        sum += delta*(x - mean);
    }    
    
    bool HasData() const
    {
        return n >= 1;
    }
    
    bool CanCalculateVariance() const
    {
        return n >= 2;
    }
    
    double Mean() const
    {
        if (HasData())
        {
            return mean;
        }
        throw "No data to calculate mean";
    }
    
    double Variance() const
    {
        if (CanCalculateVariance())
        {
            return sum / (n - 1);
        }
        throw "Not enough data to calculate variance";
    }
    
    double Deviation() const
    {
        return sqrt(Variance());
    }
    
    double Minimum() const
    {
        if (HasData())
        {
            return min;
        }
        throw "No data to calculate minimum";
    }
    
    double Maximum() const
    {
        if (HasData())
        {
            return max;
        }
        throw "No data to calculate maximum";
    }
};

//======================================================================================

const int MaxParticles = 1000;
    
void PrintUsage()
{
    std::cout <<
        "\n"
        "USAGE:\n"
        "\n"
        "fastsim fixed N dt\n"
        "    Run multiple simulations on N particles with time increment dt.\n"
        "\n";
}


int ScanNumParticles(const char *text)
{
    int n = atoi(text);
    if (n<1 || n>MaxParticles)
    {
        throw "Invalid number of particles.";
    }
    return n;
}


double ScanTimeIncrement(const char *text)
{
    double dt = atof(text);
    if (dt <= 0.0)
    {
        throw "Time increment must be a positive real number.";
    }
    return dt;
}

bool TryFixedIncrement(int n, double dt, int &maxframes)
{
    using namespace std;
    using namespace Electrons;
    
    maxframes = -1;
    
    const int NumTrials = 100;
    const int MaxFrames = 1000000;
    cout << setprecision(10) << fixed;
    StandardDeviationCalculator sd;
    for (int trial=0; trial < NumTrials; ++trial)
    {
        Simulation sim(n);
        if (sim.Converge(dt, MaxFrames))
        {
            //cout << "trial=" << trial << ", frames=" << sim.FrameNumber() << endl;
            sd.Append(static_cast<double>(sim.FrameNumber()));
        }
        else
        {
            cout << "trial=" << trial << ", *** NO CONVERGENCE ***" << endl;
            return false;
        }
    }
    
    maxframes = static_cast<int>(sd.Maximum());
    
    //cout << endl;
    cout << 
        "frames: count=" << sd.NumData() <<
        ", dt=" << setprecision(6) << dt << setprecision(2) <<
        ", mean=" << sd.Mean() << 
        ", stdev=" << sd.Deviation() << 
        ", min=" << static_cast<int>(sd.Minimum()) <<
        ", max=" << maxframes <<
        endl;        
    
    return true;
}

//======================================================================================

int AutoConverge(int n)
{
    using namespace std;
    using namespace Electrons;
    
    Simulation sim(n);
    if (sim.AutoConverge())
    {
        return 0;
    }
    return 1;
}

//======================================================================================

int main(int argc, const char *argv[])
{
    using namespace std;
    using namespace Electrons;
    try 
    {
        if (argc > 1)
        {
            const char *verb = argv[1];
            
            if (!strcmp(verb, "auto") && (argc == 3))
            {
                int n = ScanNumParticles(argv[2]);
                return AutoConverge(n);
            }
            
            if (!strcmp(verb, "fixed") && (argc == 4))
            {
                int maxframes = 0;
                int n = ScanNumParticles(argv[2]);
                double dt = ScanTimeIncrement(argv[3]);
                if (TryFixedIncrement(n, dt, maxframes))
                {
                    return 0;
                }
                return 1;
            }
            
            if (!strcmp(verb, "range") && (argc == 6))
            {
                int n = ScanNumParticles(argv[2]);
                double dt1 = ScanTimeIncrement(argv[3]);
                double dt2 = ScanTimeIncrement(argv[4]);
                double inc = ScanTimeIncrement(argv[5]);
                int bestframes = 0;
                double bestdt = 0;
                for (double dt = dt1; dt <= dt2; dt += inc)
                {
                    int maxframes = 0;
                    if (TryFixedIncrement(n, dt, maxframes))
                    {
                        if (bestframes==0 || maxframes<bestframes)
                        {
                            cout << "BEST: dt=" << setprecision(6) << dt << ", maxframes=" << maxframes << endl;
                            bestframes = maxframes;
                            bestdt = dt;
                        }
                    }
                }
                
                if (bestframes > 0)
                {
                    cout << endl;
                    cout << "Final answer: dt=" << setprecision(6) << bestdt << ", frames=" << bestframes << endl;
                }
                return 0;
            }
        }
    
        PrintUsage();
        return 2;
    }
    catch (const char *message)
    {
        cerr << "EXCEPTION: " << message << endl;
        return 3;
    }
}

