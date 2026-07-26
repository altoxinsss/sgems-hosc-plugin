#include "WeightedMomentCore.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace weighted_moment_core {
namespace {

double binomial(int n, int k)
{
    if(k < 0 || k > n)
        return 0.0;
    if(k == 0 || k == n)
        return 1.0;
    if(k > n - k)
        k = n - k;
    long double result = 1.0;
    for(int i = 1; i <= k; ++i)
    {
        result *= static_cast<long double>(n - k + i);
        result /= static_cast<long double>(i);
    }
    return static_cast<double>(result);
}

double factorial_value(int n)
{
    double value = 1.0;
    for(int i = 2; i <= n; ++i)
        value *= static_cast<double>(i);
    return value;
}

void validate_multiplicities(const int* multiplicities, int variable_count)
{
    if(multiplicities == 0 || variable_count <= 0)
        throw std::invalid_argument("weighted cumulant requires a non-empty multiplicity index.");
    for(int i = 0; i < variable_count; ++i)
    {
        if(multiplicities[i] < 0)
            throw std::invalid_argument("weighted cumulant multiplicities must be non-negative.");
    }
}

int joint_order(const int* multiplicities, int variable_count)
{
    int order = 0;
    for(int i = 0; i < variable_count; ++i)
        order += multiplicities[i];
    return order;
}

std::vector<int> labels_from_multiplicities(const int* multiplicities, int variable_count)
{
    std::vector<int> labels;
    labels.reserve(static_cast<size_t>(joint_order(multiplicities, variable_count)));
    for(int variable = 0; variable < variable_count; ++variable)
    {
        for(int repeat = 0; repeat < multiplicities[variable]; ++repeat)
            labels.push_back(variable);
    }
    return labels;
}

std::vector<double> weighted_variable_centers(
    const WeightedSample* samples,
    int sample_count,
    int variable_count,
    const WeightedMomentOptions& options)
{
    std::vector<double> centers(static_cast<size_t>(variable_count), 0.0);
    std::vector<long double> sums(static_cast<size_t>(variable_count), 0.0);
    long double weight_sum = 0.0;
    for(int i = 0; i < sample_count; ++i)
    {
        const double w = compute_weight(samples[i], options);
        if(w <= 0.0)
            continue;
        if(samples[i].values == 0 || samples[i].value_count < variable_count)
            throw std::invalid_argument("weighted cumulant sample does not contain all requested variables.");
        weight_sum += static_cast<long double>(w);
        for(int variable = 0; variable < variable_count; ++variable)
            sums[static_cast<size_t>(variable)] += static_cast<long double>(w) * samples[i].values[variable];
    }
    if(weight_sum <= 0.0)
        return centers;
    for(int variable = 0; variable < variable_count; ++variable)
        centers[static_cast<size_t>(variable)] = static_cast<double>(sums[static_cast<size_t>(variable)] / weight_sum);
    return centers;
}

double product_for_mask(
    const WeightedSample& sample,
    const std::vector<int>& labels,
    int mask,
    const std::vector<double>& centers)
{
    long double product = 1.0;
    for(int i = 0; i < static_cast<int>(labels.size()); ++i)
    {
        if((mask & (1 << i)) == 0)
            continue;
        const int variable = labels[static_cast<size_t>(i)];
        if(sample.values == 0 || variable < 0 || variable >= sample.value_count)
            throw std::invalid_argument("weighted cumulant sample does not contain all requested variables.");
        const double center = centers.empty() ? 0.0 : centers[static_cast<size_t>(variable)];
        product *= static_cast<long double>(sample.values[variable] - center);
    }
    return static_cast<double>(product);
}

void enumerate_partitions_rec(
    int n,
    int next,
    std::vector<int>& block_for_item,
    int block_count,
    std::vector< std::vector<int> >& partitions)
{
    if(next == n)
    {
        partitions.push_back(block_for_item);
        return;
    }
    for(int b = 0; b < block_count; ++b)
    {
        block_for_item[static_cast<size_t>(next)] = b;
        enumerate_partitions_rec(n, next + 1, block_for_item, block_count, partitions);
    }
    block_for_item[static_cast<size_t>(next)] = block_count;
    enumerate_partitions_rec(n, next + 1, block_for_item, block_count + 1, partitions);
}

std::vector< std::vector<int> > enumerate_partitions(int n)
{
    std::vector< std::vector<int> > partitions;
    if(n <= 0)
        return partitions;
    std::vector<int> block_for_item(static_cast<size_t>(n), 0);
    enumerate_partitions_rec(n, 1, block_for_item, 1, partitions);
    return partitions;
}

double weighted_block_moment(
    const WeightedSample* samples,
    int sample_count,
    const std::vector<int>& labels,
    int mask,
    const WeightedMomentOptions& options)
{
    if(mask == 0)
        return 1.0;
    int variable_count = 0;
    for(size_t i = 0; i < labels.size(); ++i)
        if(labels[i] + 1 > variable_count) variable_count = labels[i] + 1;
    const std::vector<double> centers = options.zero_first_order_moment
        ? weighted_variable_centers(samples, sample_count, variable_count, options)
        : std::vector<double>();

    long double weighted_sum = 0.0;
    long double weight_sum = 0.0;
    for(int i = 0; i < sample_count; ++i)
    {
        const double w = compute_weight(samples[i], options);
        if(w <= 0.0)
            continue;
        weighted_sum += static_cast<long double>(w) * product_for_mask(samples[i], labels, mask, centers);
        weight_sum += static_cast<long double>(w);
    }
    if(weight_sum <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(options.normalize_weights ? weighted_sum / weight_sum : weighted_sum);
}

}

double compute_weight(const WeightedSample& sample, const WeightedMomentOptions& options)
{
    double weight = 1.0;
    const double eps = std::numeric_limits<double>::epsilon();
    const double bandwidth = options.bandwidth > eps ? options.bandwidth : eps;
    if(options.mode == WeightMode_CoordinateGaussian)
    {
        const double dx = sample.x - options.center_x;
        const double dy = sample.y - options.center_y;
        const double dz = sample.z - options.center_z;
        const double d2 = dx * dx + dy * dy + dz * dz;
        weight = std::exp(-0.5 * d2 / (bandwidth * bandwidth));
    }
    else if(options.mode == WeightMode_ValueGaussian)
    {
        if(sample.values == 0 || sample.value_count <= 0)
            return 0.0;
        const double dv = sample.values[0] - options.reference_value;
        weight = std::exp(-0.5 * dv * dv / (bandwidth * bandwidth));
    }
    if(weight < options.min_weight)
        weight = 0.0;
    return weight;
}

double weighted_moment(
    const WeightedSample* samples,
    int sample_count,
    const int* multiplicities,
    int variable_count,
    const WeightedMomentOptions& options)
{
    if(samples == 0 || sample_count <= 0)
        throw std::invalid_argument("weighted moment requires at least one sample.");
    validate_multiplicities(multiplicities, variable_count);
    const std::vector<int> labels = labels_from_multiplicities(multiplicities, variable_count);
    if(labels.empty())
        return 1.0;
    const int full_mask = (1 << static_cast<int>(labels.size())) - 1;
    return weighted_block_moment(samples, sample_count, labels, full_mask, options);
}

double weighted_cumulant_from_multiplicities(
    const WeightedSample* samples,
    int sample_count,
    const int* multiplicities,
    int variable_count,
    const WeightedMomentOptions& options)
{
    if(samples == 0 || sample_count <= 0)
        throw std::invalid_argument("weighted cumulant requires at least one sample.");
    validate_multiplicities(multiplicities, variable_count);
    const std::vector<int> labels = labels_from_multiplicities(multiplicities, variable_count);
    const int order = static_cast<int>(labels.size());
    if(order == 0)
        return 0.0;
    if(order >= static_cast<int>(sizeof(int) * 8 - 1))
        throw std::invalid_argument("weighted cumulant order is too large for partition masks.");

    const std::vector< std::vector<int> > partitions = enumerate_partitions(order);
    long double cumulant = 0.0;
    for(size_t p = 0; p < partitions.size(); ++p)
    {
        const std::vector<int>& partition = partitions[p];
        int block_count = 0;
        for(size_t i = 0; i < partition.size(); ++i)
        {
            if(partition[i] + 1 > block_count)
                block_count = partition[i] + 1;
        }

        long double term = factorial_value(block_count - 1);
        if((block_count - 1) % 2 != 0)
            term = -term;

        for(int b = 0; b < block_count; ++b)
        {
            int mask = 0;
            for(int i = 0; i < order; ++i)
            {
                if(partition[static_cast<size_t>(i)] == b)
                    mask |= (1 << i);
            }
            term *= weighted_block_moment(samples, sample_count, labels, mask, options);
        }
        cumulant += term;
    }
    return static_cast<double>(cumulant);
}

double weighted_univariate_moment(
    const double* values,
    const double* weights,
    int sample_count,
    int order,
    bool normalize_weights,
    bool zero_first_order_moment)
{
    if(values == 0 || sample_count <= 0 || order < 0)
        throw std::invalid_argument("weighted univariate moment received invalid input.");
    if(order == 0)
        return 1.0;
    long double weighted_sum = 0.0;
    long double weight_sum = 0.0;
	long double weighted_mean_sum = 0.0;
	if(zero_first_order_moment)
	{
		for(int i = 0; i < sample_count; ++i)
		{
			const double w = (weights == 0) ? 1.0 : weights[i];
			if(w <= 0.0) continue;
			weighted_mean_sum += static_cast<long double>(w) * values[i];
			weight_sum += static_cast<long double>(w);
		}
		if(weight_sum <= 0.0)
			return std::numeric_limits<double>::quiet_NaN();
	}
	const double center = zero_first_order_moment
		? static_cast<double>(weighted_mean_sum / weight_sum)
		: 0.0;
	weight_sum = 0.0;
    for(int i = 0; i < sample_count; ++i)
    {
        const double w = (weights == 0) ? 1.0 : weights[i];
        if(w <= 0.0)
            continue;
        weighted_sum += static_cast<long double>(w) * std::pow(values[i] - center, order);
        weight_sum += static_cast<long double>(w);
    }
    if(weight_sum <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(normalize_weights ? weighted_sum / weight_sum : weighted_sum);
}

double weighted_univariate_cumulant(
    const double* values,
    const double* weights,
    int sample_count,
    int order,
    bool normalize_weights,
    bool zero_first_order_moment)
{
    if(order <= 0)
        return 0.0;
    std::vector<double> moments(static_cast<size_t>(order) + 1, 0.0);
    std::vector<double> cumulants(static_cast<size_t>(order) + 1, 0.0);
    for(int n = 1; n <= order; ++n)
    {
        moments[static_cast<size_t>(n)] = weighted_univariate_moment(
            values, weights, sample_count, n, normalize_weights, zero_first_order_moment);
    }
    cumulants[1] = moments[1];
    for(int n = 2; n <= order; ++n)
    {
        long double value = moments[static_cast<size_t>(n)];
        for(int j = 1; j < n; ++j)
        {
            value -= static_cast<long double>(binomial(n - 1, j - 1)) *
                static_cast<long double>(cumulants[static_cast<size_t>(j)]) *
                static_cast<long double>(moments[static_cast<size_t>(n - j)]);
        }
        cumulants[static_cast<size_t>(n)] = static_cast<double>(value);
    }
    return cumulants[static_cast<size_t>(order)];
}

}
