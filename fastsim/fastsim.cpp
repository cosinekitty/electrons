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
        int loop;
    
    public:
        Simulation(int numPoints)       // create an initial random state
            : frame(0)
            , loop(0)
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
        
        int FrameCount() const
        {
            return frame;
        }
        
        int LoopCount() const
        {
            return loop;
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
        
        bool AutoConverge(double shrink)
        {
            using namespace std;
            
            if (shrink<=0.0 || shrink>=1.0)
                throw "Shrink factor must be between 0 and 1.";
        
            // Theory: the simulation converges if total potential energy always decreases.
            
            // Create an auxiliary particle list to hold candidate next frames.
            const ParticleList::size_type n = particles.size();
            ParticleList nextlist;
            nextlist.reserve(n);
            for (ParticleList::size_type i = 0; i < n; ++i)
            {
                nextlist.push_back(Particle());
            }
            
            // Start out with a very large dt. Make it smaller until potential energy decreases.
            double dt = 4.0;    // known to optimally converge the n=2 case
            
            const double DeltaPowerTolerance = 1.0e-30;
            CalcTangentialForces(particles);
            double energy = PotentialEnergy(particles);
            while (true)    // frame loop: each iteration updates the particles' positions
            {
                ++frame;
                //cout << "frame=" << frame << setprecision(12) << ", energy=" << energy << ", dt=" << dt << endl;
                double nextenergy;
                while (true)    // dt adjustment loop: adjust dt as needed for potential energy to decrease
                {
                    ++loop;
                    UpdatePositions(particles, nextlist, dt);
                    CalcTangentialForces(nextlist);
                    nextenergy = PotentialEnergy(nextlist);
                    
                    // We want to detect convergence based on potential energy settling down.
                    // But because dt can change, this alone could cause dE to appear very small.
                    // So calculate the instantaneous power dP = dE/dt as a metric of convergence.
                    double dP = (nextenergy - energy) / dt;
                    //cout << "dP = " << scientific << setprecision(15) << dP << endl;
                    
                    if (fabs(dP) < DeltaPowerTolerance)
                    {
                        //cout << "final dt = " << dt << endl;
                        swap(particles, nextlist);  // get that last little bit of refinement
                        return true;    // the simulation has settled down enough!
                    }
                        
                    if (nextenergy < energy) 
                        break;
                    
                    dt *= shrink;
                    
#if 0                   
                    cout << "Decreased dt=" << setprecision(6) << fixed << dt << setprecision(12) <<
                        ", energy=" << energy << 
                        ", nextenergy=" << nextenergy << 
                        ", dE=" << scientific << (nextenergy - energy) <<
                        endl;
#endif
                        
                    if (dt < 1.0e-20)
                        return false;   // dt has decreased so much that we are probably stuck!
                }
                    
                swap(particles, nextlist);
                energy = nextenergy;
            }
        }
        
    private:
        static double CalcTangentialForces(ParticleList& particles)
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
            double maxForceMag = 0.0;
            for (Particle& p : particles)
            {
                // Calculate radial component using dot product and subtract
                // to get tangential component.
                p.force -= Vector::Dot(p.force, p.position) * p.position;
                double forceMag = p.force.Mag();
                if (forceMag > maxForceMag)
                {
                    maxForceMag = forceMag;
                }
            }
            
            return maxForceMag;
        }
        
        static void UpdatePositions(ParticleList& inlist, ParticleList& outlist, double dt)
        {
            const ParticleList::size_type numParticles = inlist.size();
            if (numParticles != outlist.size())
            {
                throw "Inconsistent particle list sizes.";
            }
            
            for (ParticleList::size_type i=0; i < numParticles; ++i)
            {
                Particle& p = inlist[i];
                Particle& q = outlist[i];                
                q.position = ((dt * p.force) + p.position).UnitVector();
            }
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
        
        static double PotentialEnergy(const ParticleList& particles)
        {
            // Potential energy is proportional to the sum of reciprocals of
            // distances between all unique pairs of particles.
        
            double energy = 0.0;
            const ParticleList::size_type numParticles = particles.size();
            for (ParticleList::size_type i=0; i < numParticles-1; ++i)
            {
                for (ParticleList::size_type j=i+1; j < numParticles; ++j)
                {
                    energy += 1.0 / (particles[i].position - particles[j].position).Mag();
                }
            }
            return energy;
        }
    };
}

//======================================================================================

const int MaxParticles = 1000;
    
void PrintUsage()
{
    std::cout <<
        "\n"
        "USAGE:\n"
        "\n"
        "fastsim auto N shrink\n"
        "    Simulate N particles using automatic convergence algorithm.\n"
        "    The 'shrink' value must be between 0 and 1, specifying how\n"
        "    much to reduce dt each time potential energy increases.\n"
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
            
            if (!strcmp(verb, "auto") && (argc == 4))
            {
                int n = ScanNumParticles(argv[2]);
                double shrink = atof(argv[3]);
                Simulation sim(n);
                if (sim.AutoConverge(shrink))
                {
                    //sim.JsonPrint(cout, 0);
                    cout << "Converged after " << sim.FrameCount() << " frames, " << sim.LoopCount() << " loops." << endl;
                    return 0;
                }
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

