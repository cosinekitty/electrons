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
#include <string>
#include <algorithm>

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
        Particle(const Vector& _position): position(_position.UnitVector()), force() {}
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
    
    // Distance Spectrum ---------------------------------------------------------
    
    struct Pair
    {
        int     aIndex;
        int     bIndex;
        double  distance;
        
        Pair(int _aIndex, int _bIndex, double _distance)
            : aIndex(_aIndex)
            , bIndex(_bIndex)
            , distance(_distance)
        {
        }
        
        void Print(std::ostream& output) const
        {
            using namespace std;
            output << 
                setw(5) << aIndex << 
                setw(5) << bIndex << 
                setw(12) << fixed << setprecision(5) << distance << 
                endl;
        }
    };
    
    bool operator < (const Pair& a, const Pair& b)      // needed for sorting
    {
        return a.distance < b.distance;
    }
    
    typedef std::vector<Pair> PairList;
    
    void Print(std::ostream& output, const PairList& pairs)
    {
        output << "DistanceSpectrum " << pairs.size() << std::endl;
        for (const Pair& p : pairs)
        {
            p.Print(output);
        }
    }
    
    // Simulation ---------------------------------------------------------------

    class Simulation
    {
    private:
        ParticleList particles;
        int frameCount;
        int updateCount;
        double energy;
    
    public:
        static const int MinParticles = 2;
        static const int MaxParticles = 1000;
        
        Simulation(int numPoints)       // create an initial random state
            : frameCount(0)
            , updateCount(0)
        {
            using namespace std;
            
            if ((numPoints < MinParticles) || (numPoints > MaxParticles))
            {
                throw "Invalid number of particles.";        
            }
                
            ifstream infile("/dev/urandom", ios::in | ios::binary);
            if (!infile) 
            {
                throw "Could not open /dev/urandom to obtain random numbers.";
            }
                
            particles.reserve(static_cast<ParticleList::size_type>(numPoints));
            for (int i=0; i < numPoints; ++i)
            {
                particles.push_back(Particle(RandomSpherePoint(infile)));
            }
            
            energy = CalcTangentialForces(particles);
        }
        
        Simulation(const char *inJsonFile)
            : frameCount(0)
            , updateCount(0)
            , energy(0)
        {
            JsonLoad(inJsonFile);
        }
        
        int ParticleCount() const
        {
            return static_cast<int>(particles.size());
        }
        
        int FrameCount() const
        {
            return frameCount;
        }
        
        int UpdateCount() const
        {
            return updateCount;
        }
        
        double PotentialEnergy() const
        {
            return energy;
        }
        
        void JsonPrint(std::ostream& output, int indent) const
        {
            using namespace std;
            
            JsonIndent(output, indent);
            output << setprecision(12) <<
                "{\"frame\":" << frameCount << 
                ", \"update\":" << updateCount << 
                ", \"energy\":" << energy << 
                ", \"particles\": [\n";
                
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
        
        bool Converge()
        {
            using namespace std;
            
            // Create auxiliary particle lists to hold candidate next frames.
            ParticleList nextlist = CreateParticleList();
            ParticleList bestlist = CreateParticleList();

            // The potential energy and tangential forces have already been calculated
            // for the current configuration.  See constructor.
                        
            while (true)    // frame loop: each iteration updates the particles' positions
            {                
                // Search for dt value that decreases potential energy as much as possible.
                // This is a balance between searching as few dt values as possible
                // and homing in on as good value as possible.
                // Calculate reasonable upper bound for dt.
                double dt = DeltaTimeUpperLimit();
                
                // Use dtAttempt to simulate a possible frame.
                double nextenergy = Update(particles, nextlist, dt);
                while (nextenergy >= energy)
                {
                    if (dt < 1.0e-20)
                    {
                        // We have made dt extremely small, yet we still can't decrease potential energy
                        // of the system. Consider the system to be converged.
                        return true;
                    }
                    dt *= 0.5;      // shrink dt and try adjusting the particle locations again
                    nextenergy = Update(particles, nextlist, dt);
                }
                
                swap(particles, nextlist);
                energy = nextenergy;
                ++frameCount;
            }
        }
        
        PairList Spectrum() const
        {
            // Generate the raw list of pair data (aIndex, bIndex, distance).
            PairList pairs;
            ParticleList::size_type numParticles = particles.size();
            pairs.reserve((numParticles * (numParticles-1)) / 2);
            for (ParticleList::size_type i=0; i < numParticles-1; ++i)
            {
                for (ParticleList::size_type j=i+1; j < numParticles; ++j)
                {
                    double distance = (particles[i].position - particles[j].position).Mag();
                    pairs.push_back(Pair(i, j, distance));
                }
            }

            // Sort the list of pairs in ascending order of distance.
            std::sort(pairs.begin(), pairs.end());
            
            return pairs;
        }
        
    private:
        double Update(ParticleList& currlist, ParticleList& nextlist, double dt)
        {
            ++updateCount;
            UpdatePositions(currlist, nextlist, dt);
            return CalcTangentialForces(nextlist);      // returns updated potential energy of particles in 'nextlist'
        }
    
        double DeltaTimeUpperLimit() const
        {
            // Goal: find a reasonable upper bound for simulation increment dt.
            // Strategy: Approximate the value of dt that would result in the closest
            // pair of particles (ignoring all other particles) to increase their
            // distance from each other by no more than their current distance times alpha.
            const double alpha = 0.5;
            
            // Find the minimum distance between any pair of particles.
            double r = MinimumPairDistance(particles);
            
            // dr = 2*dt/(r^2), because both particles move by dt/(r^2).
            // dr = alpha*r.
            // Therefore, dt = (alpha/2)*r^3.
            
            return (alpha/2.0) * (r*r*r);
        }
        
        static double MinimumPairDistance(const ParticleList& particles)
        {
            const ParticleList::size_type numParticles = particles.size();
            if (numParticles < 2)
            {
                throw "There must be at least 2 particles to find a minimum pair distance.";
            }
            
            double minDistanceSquared = 5.0;  // particles can never be more than 2 units apart, and 2^2 = 4.
            for (ParticleList::size_type i=0; i < numParticles-1; ++i)
            {
                for (ParticleList::size_type j=i+1; j < numParticles; ++j)
                {
                    double distanceSquared = (particles[i].position - particles[j].position).MagSquared();
                    if (distanceSquared < minDistanceSquared)
                    {
                        minDistanceSquared = distanceSquared;
                    }
                }
            }
            
            if (minDistanceSquared > 4.000000001)       // allow ample room for roundoff error in the n=2 case
            {
                throw "Impossible minimum distance between particles!";
            }
            
            return sqrt(minDistanceSquared);
        }
    
        static double CalcTangentialForces(ParticleList& particles)
        {
            // Reset each particle's force to a zero vector.
            for (Particle& p : particles)
            {
                p.force.Reset();
            }
            
            // Compute force between each unique pair of particles
            // and total potential energy of the system.
            double energy = 0.0;
            const ParticleList::size_type numParticles = particles.size();
            for (ParticleList::size_type i=0; i < numParticles-1; ++i)
            {
                for (ParticleList::size_type j=i+1; j < numParticles; ++j)
                {
                    energy += AddForces(particles[i], particles[j]);
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
            
            return energy;
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
                // IMPORTANT: q.force is left whatever random garbage was there before - NOT VALID VALUES!
            }
        }
        
        static double AddForces(Particle& a, Particle& b)
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
            
            // Calculate potential energy of this pair and return it.
            // Potential energy is proportional to 1/r = sqrt(1/r^2).
            return sqrt(forcemag);
        }
        
        ParticleList CreateParticleList()
        {
            const ParticleList::size_type n = particles.size();
            ParticleList list;
            list.reserve(n);
            for (ParticleList::size_type i = 0; i < n; ++i)
            {
                list.push_back(Particle());
            }
            return list;
        }
        
        static const int MAXSTRING = 40;

        void AppendChar(char *string, int& index, char c)
        {
            if (index < MAXSTRING)
            {
                string[index++] = c;
                string[index] = '\0';
            }
            else
            {
                throw "JSON token is too long.";
            }
        }

        void JsonLoad(const char *inFileName)
        {
            using namespace std;
            
            // WARNING: This is not a real JSON parser. It is a hack that takes advantage
            // of the format we know we write to JSON files.
            // If you need this program to load JSON that you create, try to follow the format
            // that this program outputs as closely as possible.
            // Strategy: we only load the data we can't recalculate (counters and position vectors).
            // We recalculate energy and force vectors using the position vectors.
            
            ifstream infile(inFileName);
            if (!infile)
            {
                throw "Cannot open input file.";
            }
            
            particles.clear();
            frameCount = 0;
            updateCount = 0;
            
            // Here is the hack: always remember the last 2 strings surrounded by double quotes.
            // Then look for all sequences of digits, minus signs, or periods.
            // We can determine what each number means by looking at the last 1 or 2 strings.

            char key[1 + MAXSTRING];
            char numeric[1 + MAXSTRING];
            key[0] = numeric[0] = '\0';
            int k = 0;   // index into key
            int n = 0;   // index into numeric
            char c;
            enum {SEARCH, KEY, NUMBER} state = SEARCH;
            bool inPositionVector = false;
            Vector vec;

            while ((infile >> c).good())
            {
                switch (state)
                {
                    case SEARCH:
                        if ((c >= '0' && c <= '9') || c=='-')
                        {
                            state = NUMBER;
                            n = 0;
                            AppendChar(numeric, n, c);
                        }
                        else if (c == '"')
                        {
                            state = KEY;
                            k = 0;
                        }
                        break;

                    case KEY:
                        if (c == '"')
                        {
                            state = SEARCH;
                            if (!strcmp(key, "position"))
                            {
                                inPositionVector = true;
                            }
                            else if (key[0]<'x' || key[0]>'z' || key[1])
                            {
                                if (inPositionVector)
                                {
                                    inPositionVector = false;
                                    particles.push_back(Particle(vec));
                                    vec.Reset();
                                }
                            }
                        }
                        else
                        {
                            AppendChar(key, k, c);
                        }
                        break;

                    case NUMBER:
                        if ((c >= '0' && c <= '9') || c=='-' || c=='+' || c=='e' || c=='E' || c=='.')
                        {
                            AppendChar(numeric, n, c);
                        }
                        else
                        {
                            state = SEARCH;
                            double value = atof(numeric);
                            if (inPositionVector)
                            {
                                if (key[0] && !key[1])
                                {
                                    switch (key[0])
                                    {
                                        case 'x': vec.x = value; break;
                                        case 'y': vec.y = value; break;
                                        case 'z': vec.z = value; break;
                                    }
                                }
                            }
                            else if (!strcmp(key, "frame"))
                            {
                                frameCount = atoi(numeric);
                            }
                            else if (!strcmp(key, "update") || !strcmp(key, "loop"))
                            {
                                updateCount = atoi(numeric);
                            }
                        }
                        break;

                    default:
                        throw "Unknown state trying to load JSON";
                }
            } 
            
            energy = CalcTangentialForces(particles);
        }
    };
}

//======================================================================================
    
void PrintUsage()
{
    std::cout <<
        "\n"
        "USAGE: (where N=" << Electrons::Simulation::MinParticles << ".." << Electrons::Simulation::MaxParticles << ")\n"
        "\n"
        "fastsim converge N [outfile.json]\n"
        "    Simulate N particles until they reach minimum potential energy.\n"
        "    If output filename is omitted, it defaults to 'sim.json'.\n"
        "\n"
        "fastsim random N outfile.json\n"
        "    Save a random configuration of N particles to the specified file.\n"
        "\n"
        "fastsim copy infile.json outfile.json\n"
        "    Loads simulation from infile, then writes to outfile.\n"
        "    Not much practical value; used as a unit test of the loader.\n"
        "\n"
        "fastsim compare a.json b.json\n"
        "    Compares the two simulations to see if they have the same\n"
        "    relative configuration of particle positions, within tolerance.\n"
        "\n"
        "fastsim spectrum infile.json\n"
        "    Calculates the distance spectrum of the given simulation,\n"
        "    defined as the list of all distances between pairs of\n"
        "    particles, sorted in ascending order.\n"
        "\n";
}

int ScanNumParticles(const char *text)
{
    int n = atoi(text);
    if ((n < Electrons::Simulation::MinParticles) || (n > Electrons::Simulation::MaxParticles))
    {
        throw "Invalid number of particles.";
    }
    return n;
}

void Save(Electrons::Simulation& sim, const char *outFileName)
{
    std::ofstream outfile(outFileName);
    if (!outfile)
    {
        throw "Error opening output file!";
    }
    sim.JsonPrint(outfile, 0);
    if (!outfile)
    {
        throw "Error writing to output file!";
    }
}

//======================================================================================

inline bool Different(double x, double y, const char *label, double tolerance)
{
    using namespace std;

    double diff = fabs(x - y);
    if (diff >= tolerance)
    {
        cout << setprecision(12) <<
            "Different " << label << 
            ": x=" << x << 
            ", y=" << y << 
            ", diff=" << diff << 
            ", tolerance=" << tolerance << endl;
            
        return true;
    }
    return false;
}

bool Compare(Electrons::Simulation& asim, Electrons::Simulation& bsim)
{
    using namespace std;

    const int n = asim.ParticleCount();
    if (n != bsim.ParticleCount())
    {
        cout << "Simulations have different particle counts." << endl;
        return false;
    }
    
    const double tolerance = 1.0e-7;
    const int pairs = (n * (n-1)) / 2;

    if (Different(asim.PotentialEnergy(), bsim.PotentialEnergy(), "potential energies", pairs*tolerance))
    {
        return false;
    }

    return true;
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
            
            if (!strcmp(verb, "converge") && (argc >= 3) && (argc <= 4))
            {
                int n = ScanNumParticles(argv[2]);
                const char *outFileName = (argc > 3) ? argv[3] : "sim.json";
                if (0 == remove(outFileName))
                {
                    cout << "Deleted existing file '" << outFileName << "'" << endl;
                }
                Simulation sim(n);
                if (sim.Converge())
                {
                    cout << "Converged after " << sim.FrameCount() << " frames, " << sim.UpdateCount() << " updates." << endl;
                    Save(sim, outFileName);
                    return 0;
                }
                cout << "FAILED TO CONVERGE!" << endl;
                return 1;
            }
            
            if (!strcmp(verb, "random") && (argc == 4))
            {
                int n = ScanNumParticles(argv[2]);
                const char *outFileName = argv[3];
                Simulation sim(n);
                Save(sim, outFileName);
                return 0;
            }

            if (!strcmp(verb, "copy") && (argc == 4))
            {
                // Test the JSON loader.
                const char *inFileName = argv[2];
                const char *outFileName = argv[3];
                Simulation sim(inFileName);
                Save(sim, outFileName);
                return 0;
            }

            if (!strcmp(verb, "compare") && (argc == 4))
            {
                const char *aFileName = argv[2];
                const char *bFileName = argv[3];
                Simulation asim(aFileName);
                Simulation bsim(bFileName);
                if (Compare(asim, bsim))
                {
                    cout << "The simulations are equivalent." << endl;
                    return 0;
                }
                return 9;   // special return value that scripts can use to find a surprising convergence!
            }
            
            if (!strcmp(verb, "spectrum") && (argc == 3))
            {
                const char *inFileName = argv[2];
                Simulation sim(inFileName);
                PairList spectrum = sim.Spectrum();
                Print(cout, spectrum);
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

