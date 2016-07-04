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


