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

#include "lodepng.h"

namespace Electrons
{
    class Simulation;
}

void Save(Electrons::Simulation& sim, const char *outFileName);

namespace Electrons
{

    // lodepng helpers...
    const unsigned BYTES_PER_PIXEL = 4;

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
        bool   visited;     // a utility flag

        Particle(): position(), force(), visited(false) {}
        Particle(const Vector& _position): position(_position.UnitVector()), force(), visited(false) {}
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

    inline std::ostream& operator << (std::ostream& output, const Vector& v)
    {
        JsonPrint(output, v);
        return output;
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
        int     group;

        Pair(int _aIndex, int _bIndex, double _distance)
            : aIndex(_aIndex)
            , bIndex(_bIndex)
            , distance(_distance)
            , group(0)
            {}

        Pair(const Pair& other, int _group)
            : aIndex(other.aIndex)
            , bIndex(other.bIndex)
            , distance(other.distance)
            , group(_group)
            {}

        void Print(std::ostream& output) const
        {
            using namespace std;
            output <<
                setw(5) << aIndex <<
                setw(5) << bIndex <<
                setw(14) << fixed << setprecision(10) << distance <<
                setw(5) << group <<
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

    struct GroupPattern
    {
        PairList pattern;
        std::vector<int> groups;    // groups[particleIndex] = groupNumber

        GroupPattern(const PairList& spectrum, double tolerance)
        {
            // Find total number of particles in the spectrum by scanning the particle indices in it.
            int n = 0;
            for (const Pair& p : spectrum)
            {
                if (p.aIndex >= n)
                {
                    n = 1 + p.aIndex;
                }
                if (p.bIndex >= n)
                {
                    n = 1 + p.bIndex;
                }
            }

            // Initialize each particle's group to 0, meaning "no group".
            groups.reserve(static_cast<std::vector<int>::size_type>(n));
            for (int i=0; i < n; ++i)
            {
                groups.push_back(0);
            }

            // Break the spectrum into bands of approximately equal pair distances.
            // Once a particle has been assigned to a group, exclude it from any later groups.
            // This is the same idea used for drawing the pattern lines in the web version.

            int g = 1;
            int groupCount = 0;
            double prevDistance = -1.0;
            for (const Pair& p : spectrum)
            {
                if (p.distance - prevDistance > tolerance)
                {
                    if (groupCount > 0)     // don't "waste" group numbers on empty groups
                    {
                        ++g;                // starting a new group number
                        groupCount = 0;     // we have not found any pairs in this group yet
                    }
                    prevDistance = p.distance;
                }

                if ((groups[p.aIndex]==0 || groups[p.aIndex]==g) && (groups[p.bIndex]==0 || groups[p.bIndex]==g))
                {
                    groups[p.aIndex] = groups[p.bIndex] = g;
                    pattern.push_back(Pair(p, g));
                    ++groupCount;
                }
            }
        }

        void Print(std::ostream& output)
        {
            using namespace std;

            Electrons::Print(output, pattern);

            output << "Groups " << groups.size() << endl;
            int i = 0;
            for (int g : groups)
            {
                output << setw(5) << i << setw(5) << g << endl;
                ++i;
            }
        }

        static bool Compatible(const PairList& a, const PairList& b, double tolerance)
        {
            int n = static_cast<int>(a.size());
            if (static_cast<int>(b.size()) != n)
                return false;

            for (int i=0; i < n; ++i)
            {
                double error = fabs(a[i].distance - b[i].distance);
                if (error > tolerance)
                    return false;
            }

            return true;
        }
    };

    // Simulation ---------------------------------------------------------------

    class Simulation
    {
    private:
        ParticleList particles;
        int frameCount;
        int updateCount;
        int votes;          // used for statistical purposes
        double energy;

    public:
        static const int MinParticles = 2;
        static const int MaxParticles = 1000;

        Simulation(int numPoints)       // create an initial random state
            : frameCount(0)
            , updateCount(0)
            , votes(0)
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
            , votes(0)
            , energy(0)
        {
            JsonLoad(inJsonFile);
        }

        void AddVote()
        {
            ++votes;
        }

        int Votes() const
        {
            return votes;
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

        void Converge()
        {
            using namespace std;

            // Create auxiliary particle lists to hold candidate next frames.
            ParticleList nextlist = CreateParticleList();

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
                        return;
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

        void Draw(const char *outFileName) const
        {
            using namespace std;

            const int pixelsWide = 500;
            const int pixelsHigh = 500;
            //const unsigned char OPAQUE_ALPHA_VALUE = 255;
            const unsigned RGBA_BUFFER_SIZE = pixelsWide * pixelsHigh * BYTES_PER_PIXEL;

            vector<unsigned char> rgbaBuffer(RGBA_BUFFER_SIZE);
            fill(rgbaBuffer.begin(), rgbaBuffer.end(), 255);    // fill image with opaque white pixels

            Render(rgbaBuffer, pixelsWide, pixelsHigh, pixelsWide/2, pixelsHigh/2, pixelsWide/2,
                Vector(1, 0, 0),
                Vector(0, 1, 0));

            unsigned error = lodepng::encode(
                outFileName,
                rgbaBuffer,
                pixelsWide,
                pixelsHigh);

            if (error != 0)
            {
                cerr << "lodepng::encode() returned error " << error << endl;
                throw "Image encoding error";
            }
        }

        void Fix(int northPoleIndex, int eastLineIndex)
        {
            double oldEnergy = energy;

            // InternalFix does not update force vectors and potential energy.
            // This makes Compare() much more efficient.
            // But doing so is required to preserve the integrity for public callers.
            InternalFix(northPoleIndex, eastLineIndex);

            UpdateAfterRotation(oldEnergy);
        }

        // The Compare() function has Simulation objects passed by value, not reference.
        // This is intentional, because the algorithm modifies the state of both asim and bsim.
        // This prevents unintentional side-effects to the caller.
        static bool Compare(Simulation asim, Simulation bsim)
        {
            using namespace std;

            const int n = asim.ParticleCount();
            if (n != bsim.ParticleCount())
            {
                throw "Simulations have different particle counts.";
            }

            if (n < 2)
            {
                // Trivial case: any pair of sims both with 0 or 1 particles
                // have identical configurations, by definition.
                // This short-cut is needed to avoid exceptions below.
                return true;
            }

            // Obtain the distance spectra of both simulations.
            const double ErrorTolerance = 1.0e-5;
            const double GroupTolerance = 1.0e-3;
            GroupPattern ag(asim.Spectrum(), GroupTolerance);
            GroupPattern bg(bsim.Spectrum(), GroupTolerance);

            // See if the two simulations have compatible spectrum groups.
            if (GroupPattern::Compatible(ag.pattern, bg.pattern, GroupTolerance))
            {
                // Use asim spectrum to find the two particles in asim that are closest together.
                // (This is not a unique answer; there may be many pairs of the same minimal length.)
                // Fix the asim sphere based on those two particles.
                int a1 = ag.pattern[0].aIndex;
                int a2 = ag.pattern[0].bIndex;
                asim.InternalFix(a1, a2);

                // Find every distinct pair of particles in bsim that are the same distance apart.
                // We already tested compatibility between the respective groups, so try every
                // pair in bsim's first group, and try each pair in both directions.
                const int numPairs = static_cast<int>(bg.pattern.size());
                for (int i=0; (i < numPairs) && (bg.pattern[i].group == bg.pattern[0].group); ++i)
                {
                    // Fix the bsim sphere based on those two particles.
                    // Calculate mean position error between asim and bsim.
                    // If the error is below tolerance, declare the simulations equivalent.
                    int b1 = bg.pattern[i].aIndex;
                    int b2 = bg.pattern[i].bIndex;

                    if (TryOrientation(asim, bsim, b1, b2) < ErrorTolerance)
                        return true;

                    if (TryOrientation(asim, bsim, b2, b1) < ErrorTolerance)
                        return true;
                }
            }
            return false;
        }

    private:
        void Mirror()
        {
            // Tricky: this mirror algorithm has to be compatible with InternalFix().
            // After calling InternalFix, the simulation leaves the first particle
            // at (0, 0, 1).  Mirror() cannot change the location of that particle.
            // The second particle is at (x, 0, z), where x >= 0.
            // That particle also must remain in the same place.
            // Both particles have y=0, so by negating all the y coordinates,
            // we create a mirror image without moving the two important points at all.

            for (Particle& p : particles)
            {
                p.position.y *= -1.0;
            }
        }

        static double TryOrientation(
            Simulation& asim,
            Simulation& bsim,
            int b1,
            int b2)
        {
            bsim.InternalFix(b1, b2);
            double error1 = MeanPositionError(asim, bsim);
            bsim.Mirror();
            double error2 = MeanPositionError(asim, bsim);
            bsim.Mirror();
            return (error1 < error2) ? error1 : error2;
        }

        void UnvisitAllParticles()
        {
            for (Particle& p : particles)
            {
                p.visited = false;
            }
        }

        static double MeanPositionError(Simulation& asim, Simulation& bsim)
        {
            if (asim.ParticleCount() != bsim.ParticleCount())
            {
                throw "Simulations must have same particle count.";
            }

            if (asim.ParticleCount() < 1)
            {
                throw "Simulations must have at least one particle each.";
            }

            // Find the root-mean-squared error of all vector positions.
            // The two simulations can have particles in any order, so we compare
            // each with the others.
            double sum = 0.0;
            bsim.UnvisitAllParticles();
            for (Particle &ap : asim.particles)
            {
                // Find the particle bp that is closest to ap.
                // Exclude bp if we have already paired it with another ap.
                Particle *closest = nullptr;
                double bestError = 5.0;      // largest possible distance is 2.0, so largest error is 4.0.

                for (Particle &bp : bsim.particles)
                {
                    if (!bp.visited)
                    {
                        double error = (ap.position - bp.position).MagSquared();
                        if (error < bestError)
                        {
                            bestError = error;
                            closest = &bp;
                        }
                    }
                }

                if (closest == nullptr)
                {
                    throw "Could not find partner particle.";
                }

                closest->visited = true;    // prevent us from using this particle from bsim again
                sum += bestError;
            }

            return sqrt(sum) / asim.ParticleCount();
        }

        void UpdateAfterRotation(double oldEnergy)
        {
            // Recalculate force vectors and potential energy for the new configuration.
            // This is so we maintain the invariant of this class, but also serves as a sanity check.
            energy = CalcTangentialForces(particles);
            if (fabs(energy - oldEnergy) > 1.0e-6)
            {
                throw "Potential energy should not have changed.";
            }
        }

        void InternalFix(int northPoleIndex, int eastLineIndex)
        {
            //
            // Rotate the sphere so that the two particles specified by index
            // are oriented in a "standard" way:
            //
            //      particle[northPoleIndex] will end up at (0, 0, 1).
            //      particle[eastLineIndex]  will end up at (x, 0, z) where x >= 0.
            //
            // This method allows us to more easily compare two Simulations
            // to see if they have the same relative configuration.
            //

            // Validate the particle indices.
            const int n = ParticleCount();
            if (n < 2)
            {
                throw "Simulation must have at least 2 particles.";
            }

            if (northPoleIndex<0 || northPoleIndex>=n)
            {
                throw "Invalid particle index to move to north pole.";
            }

            if (eastLineIndex<0 || eastLineIndex>=n)
            {
                throw "Invalid particle index to move to east line.";
            }

            if (northPoleIndex == eastLineIndex)
            {
                throw "Particle indices must be different.";
            }

            const double tolerance = 1.0e-6;

            Vector& north = particles[northPoleIndex].position;
            Vector& east = particles[eastLineIndex].position;

            // One final sanity check: the position vectors must not lie on parallel lines
            // if there are more than 2 particles.  Otherwise, there are an infinite
            // number of ways to rotate the sphere, since they will both end up on the
            // north pole, or one will end up on the north pole and the other on the south pole.
            if (n > 2)
            {
                // If the two particles have position vectors that point in the same
                // direction, or in exactly opposite directions, then their
                // dot product will be (close to) +1 or -1, respectively.
                double alignment = fabs(1.0 - fabs(Vector::Dot(north, east)));
                if (alignment < tolerance)
                {
                    // Dump coordinates for diagnostics.
                    std::cout << "north=" << north << std::endl;
                    std::cout << "east =" << east  << std::endl;
                    throw "Particles not allowed to lie on the same diameter line.";
                }
            }

            // Rotate the sphere around the z-axis to make north.y = 0 and north.x >= 0.
            RotateZ(+north.x, -north.y);

            // Now rotate around the y-axis so that north.x = 0 and north.y = 0.
            RotateY(+north.z, -north.x);

            // Rotate again around the z-axis so that the east point has y=0, x>=0.
            RotateZ(+east.x, -east.y);

            // Sanity-check that we met our goals.
            if (fabs(north.x) > tolerance || fabs(north.y) > tolerance || fabs(east.y) > tolerance)
            {
                // Dump coordinates for diagnostics.
                std::cout << "north=" << north << std::endl;
                std::cout << "east =" << east  << std::endl;
                throw "Vectors did not end up in correct locations.";
            }
        }

        void RotateY(double a, double b)
        {
            if (NormalizeRotation(a, b))
            {
                for (Particle& p : particles)
                {
                    // [See explanation in RotateZ function.]
                    // (a + ib)*(z + ix) = (a*z - b*x) + i(a*x + b*z)
                    double z = a*p.position.z - b*p.position.x;
                    p.position.x = a*p.position.x + b*p.position.z;
                    p.position.z = z;
                }
            }
        }

        void RotateZ(double a, double b)
        {
            if (NormalizeRotation(a, b))
            {
                // We rotate around the z-axis by treating each point (x, y, z)
                // as a complex number (x + iy), while leaving z constant.
                // Then multiply (x + iy) by another complex number (a + ib) where
                // a = cos(theta), y = sin(theta), theta = rotation angle.
                // We avoid needing to use cos and sin function calls by noticing that
                // cos(-theta) = x/mag, sin(-theta) = y/mag.
                // The negative in (-theta) is because we want to rotate opposite to
                // the (x, y) coordinates of the vector so that its x becomes 0.
                // In other words, by multiplying (x + iy) by its complex conjugate (x - iy)
                // dividing by its magnitude |x + iy|, we end up with |x + 0i|.
                // All other points are treated the same way and end up rotated an equivalent amount.
                for (Particle& p : particles)
                {
                    // (a + ib)*(x + iy) = (a*x - b*y) + i(a*y + b*x)
                    double x = a*p.position.x - b*p.position.y;
                    p.position.y = a*p.position.y + b*p.position.x;
                    p.position.x = x;
                }
            }
        }

        static bool NormalizeRotation(double& a, double& b)
        {
            double mag = sqrt(a*a + b*b);
            if (mag > 1.0e-6)
            {
                a /= mag;
                b /= mag;
                return true;
            }
            return false;
        }

        void Render(
            std::vector<unsigned char>& rgbaBuffer,
            int pixelsWide,
            int pixelsHigh,
            int horCenter,
            int verCenter,
            int pixelRadius,
            Vector horVector,
            Vector verVector) const
        {
            const int dotradius = 2;

            for (const Particle& p : particles)
            {
                int hh = horCenter + static_cast<int>(pixelRadius * Vector::Dot(horVector, p.position));
                int vv = verCenter + static_cast<int>(pixelRadius * Vector::Dot(verVector, p.position));
                for (int dh = -dotradius; dh <= +dotradius; ++dh)
                {
                    for (int dv = -dotradius; dv <= +dotradius; ++dv)
                    {
                        if (dh*dh + dv*dv <= dotradius*dotradius)
                        {
                            int h = hh + dh;
                            int v = vv + dv;

                            if ((h >= 0) && (h < pixelsWide) && (v >= 0) && (v <= pixelsHigh))
                            {
                                int index = BYTES_PER_PIXEL * ((v*pixelsWide) + h);
                                rgbaBuffer[index++] = 0;    // red
                                rgbaBuffer[index++] = 0;    // green
                                rgbaBuffer[index++] = 0;    // blue
                            }
                        }
                    }
                }
            }
        }

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

        ParticleList CreateParticleList() const
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

        static void AppendChar(char *string, int& index, char c)
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
        "fastsim spectrum infile.json tolerance\n"
        "    Calculates the distance spectrum of the given simulation,\n"
        "    defined as the list of all distances between pairs of\n"
        "    particles, sorted in ascending order.\n"
        "    Tolerance is a positive real number (e.g. 0.001) specifying\n"
        "    how much distance error to tolerate when grouping the particles.\n"
        "\n"
        "fastsim draw infile.json outfile.png\n"
        "    Make a png image of the given simulation.\n"
        "\n"
        "fastsim fix infile.json outfile.json north east\n"
        "    Normalize the orientation of the simulation in infile.json\n"
        "    such that the particle at index 'north' ends up at the\n"
        "    north pole (0, 0, 1) and that for 'east' ends up on an\n"
        "    easterly meridian (x >= 0, 0, z).\n"
        "\n"
        "fastsim search N1 N2 limit\n"
        "    Repeatedly generate N-particle simulations and converge them,\n"
        "    where N is increased in the closed range [N1..N2].\n"
        "    Search for all different ways each number N of particles can settle.\n"
        "    Stop searching after 'limit' attempts for each value of N.\n"
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

void Search(int numParticles, int limit)
{
    using namespace std;
    using namespace Electrons;

    vector<Simulation> simlist;         // list of distinct settling patterns we have found
    int count;

    for (count=0; count < limit; ++count)
    {
        Simulation bsim(numParticles);
        bsim.Converge();

        // Compare this settled pattern with every distinct one we know so far.
        bool found = false;
        for (Simulation& asim : simlist)
        {
            if (Simulation::Compare(asim, bsim))
            {
                asim.AddVote();
                found = true;
                break;
            }
        }

        if (!found)
        {
            bsim.AddVote();
            simlist.push_back(bsim);
            string outFileName = "search" + to_string(numParticles) + "-" + to_string(simlist.size()) + ".json";
            Save(bsim, outFileName.c_str());
            cout << "Saved " << outFileName << " : energy=" << bsim.PotentialEnergy() << ", loop=" << count << endl;
        }
    }

    cout << endl << "Done!  Vote tallies:" << endl;
    count = 0;
    for (Simulation &sim : simlist)
    {
        ++count;
        cout << sim.Votes() << " votes for #" << count << endl;
    }
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
                sim.Converge();
                cout << "Converged after " << sim.FrameCount() << " frames, " << sim.UpdateCount() << " updates." << endl;
                Save(sim, outFileName);
                return 0;
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
                if (Simulation::Compare(asim, bsim))
                {
                    cout << "The simulations are equivalent." << endl;
                    return 0;
                }
                cout << "SIMULATIONS DID NOT MATCH!" << endl;
                return 9;   // special return value that scripts can use to find a surprising convergence!
            }

            if (!strcmp(verb, "spectrum") && (argc == 4))
            {
                const char *inFileName = argv[2];
                double tolerance = atof(argv[3]);
                if (tolerance <= 0.0)
                {
                    throw "Tolerance must be a positive real number.";
                }
                Simulation sim(inFileName);
                PairList spectrum = sim.Spectrum();
                //Print(cout, spectrum);
                GroupPattern gp(spectrum, tolerance);
                gp.Print(cout);
                return 0;
            }

            if (!strcmp(verb, "draw") && (argc == 4))
            {
                const char *inFileName = argv[2];
                const char *outFileName = argv[3];
                Simulation sim(inFileName);
                sim.Draw(outFileName);
                return 0;
            }

            if (!strcmp(verb, "fix") && (argc == 6))
            {
                const char *inFileName = argv[2];
                const char *outFileName = argv[3];
                int northPoleIndex = atoi(argv[4]);
                int eastLineIndex = atoi(argv[5]);
                Simulation sim(inFileName);
                sim.Fix(northPoleIndex, eastLineIndex);
                Save(sim, outFileName);
                return 0;
            }

            if (!strcmp(verb, "search") && (argc == 5))
            {
                int n1 = ScanNumParticles(argv[2]);
                int n2 = ScanNumParticles(argv[3]);
                int limit = atoi(argv[4]);
                for (int n = n1; n <= n2; ++n)
                {
                    Search(n, limit);
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
