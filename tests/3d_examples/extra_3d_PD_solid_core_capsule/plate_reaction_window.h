#ifndef PLATE_REACTION_WINDOW_H
#define PLATE_REACTION_WINDOW_H

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>

class PlateReactionWindow
{
  public:
    explicit PlateReactionWindow(double duration) : duration_(duration) {}

    void observe(double time, double upper, double lower)
    {
        samples_.push_back({time, upper, lower});
        const double cutoff = time - duration_;
        while (samples_.size() > 1 && samples_[1].time <= cutoff)
            samples_.pop_front();
    }

    double relativeSpan() const
    {
        if (samples_.empty() || samples_.front().time > samples_.back().time - duration_)
            return std::numeric_limits<double>::infinity();

        double low = std::numeric_limits<double>::infinity();
        double high = 0;
        for (const auto &sample : samples_)
        {
            if (!std::isfinite(sample.upper) || !std::isfinite(sample.lower) ||
                sample.upper <= 0 || sample.lower <= 0)
                return std::numeric_limits<double>::infinity();
            low = std::min({low, sample.upper, sample.lower});
            high = std::max({high, sample.upper, sample.lower});
        }
        return (high - low) / (0.5 * (high + low));
    }

    void clear() { samples_.clear(); }

  private:
    struct Sample
    {
        double time;
        double upper;
        double lower;
    };

    double duration_;
    std::deque<Sample> samples_;
};

#endif
