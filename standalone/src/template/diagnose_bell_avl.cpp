#include "CumulantCal.h"
#include "Initialization.h"
#include "SearchNode.h"
#include "config.h"
#include "grid_io.h"
#include "text.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static std::string get_alias(
    const without_sgems::IniConfig& config,
    const std::string& section,
    const std::vector<std::string>& keys)
{
    for(const std::string& key : keys)
    {
        if(config.has(section, key))
        {
            return config.get(section, key);
        }
    }
    throw std::runtime_error("missing config value");
}

static std::string get_alias_or(
    const without_sgems::IniConfig& config,
    const std::string& section,
    const std::vector<std::string>& keys,
    const std::string& fallback)
{
    for(const std::string& key : keys)
    {
        if(config.has(section, key))
        {
            return config.get(section, key);
        }
    }
    return fallback;
}

static fs::path resolve_path(const fs::path& config_path, const std::string& value)
{
    fs::path path(value);
    if(path.is_absolute())
    {
        return path;
    }
    return config_path.parent_path() / path;
}

static void free_keys(Keys* key)
{
    if(key == NULL)
    {
        return;
    }
    std::free(key->NOD);
    std::free(key->NLAG);
    std::free(key->LAG);
    std::free(key->LTOL);
    std::free(key->AZM);
    std::free(key->ATOL);
    std::free(key->BANDWH);
    std::free(key->DIP);
    std::free(key->DTOL);
    std::free(key->BANDWD);
    std::free(key);
}

static Keys* build_keys(
    const without_sgems::GridData& grid,
    const std::vector<int>& nod,
    const std::vector<int>& nlag,
    const std::vector<int>& lag,
    const std::vector<int>& ltol,
    const std::vector<int>& azm,
    const std::vector<int>& atol,
    const std::vector<int>& bandwh,
    const std::vector<int>& dip,
    const std::vector<int>& dtol,
    const std::vector<int>& bandwd)
{
    const size_t dir_count = nod.size();
    Keys* key = static_cast<Keys*>(std::calloc(1, sizeof(Keys)));
    key->DIM = (grid.nz <= 1) ? 2 : 3;
    key->REG_IRREG = 1;
    key->NX = grid.nx;
    key->NY = grid.ny;
    key->NZ = grid.nz;
    key->DIMX = (grid.nx > 0) ? grid.nx - 1 : 0;
    key->DIMY = (grid.ny > 0) ? grid.ny - 1 : 0;
    key->DIMZ = (grid.nz > 0) ? grid.nz - 1 : 0;
    key->ORDER = static_cast<int>(dir_count) + 1;
    key->NDIR = static_cast<int>(dir_count);
    key->NOD = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->NLAG = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->LAG = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->LTOL = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->AZM = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->ATOL = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->BANDWH = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->DIP = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->DTOL = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->BANDWD = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    for(size_t i = 0; i < dir_count; ++i)
    {
        key->NOD[i] = nod[i];
        key->NLAG[i] = nlag[i];
        key->LAG[i] = lag[i];
        key->LTOL[i] = ltol[i];
        key->AZM[i] = azm[i];
        key->ATOL[i] = atol[i];
        key->BANDWH[i] = bandwh[i];
        key->DIP[i] = dip[i];
        key->DTOL[i] = dtol[i];
        key->BANDWD[i] = bandwd[i];
    }
    return key;
}

static int make_index(int full_order, const std::vector<int>& labels, std::vector<int>& out)
{
    out.assign(static_cast<size_t>(full_order) + 1, 0);
    out[0] = static_cast<int>(labels.size());
    for(int label : labels)
    {
        if(label < 1 || label > full_order)
        {
            return -1;
        }
        out[static_cast<size_t>(label)]++;
    }
    return 0;
}

static double avl_eval(int* index, int order)
{
    SetFactorSwitchesForRun(1, 1, 1, 1);
    if(PrecomputeMomentsForIndex(index, order) == myEmpty)
    {
        return myEmpty;
    }
    if(InitRoot2(order) == -1)
    {
        return myEmpty;
    }
    return BetterRecCumulant3(index, order);
}

int main(int argc, char** argv)
{
    try
    {
        fs::path config_path = argc > 1 ? fs::path(argv[1]) : fs::path("z25_hosc_avl_default.ini");
        int lag_index = argc > 2 ? std::atoi(argv[2]) : 17;
        const without_sgems::IniConfig config = without_sgems::IniConfig::load(config_path);
        const fs::path source_grid_path = resolve_path(
            config_path, get_alias(config, "general", {"source_grid", "source_grid_", "Source_Grid"}));
        const without_sgems::GridData grid = without_sgems::load_hosc_grid(source_grid_path);

        Keys* key = build_keys(
            grid,
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"number_of_occurrences_in_each_direction", "nod"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"number_of_lags_in_each_direction", "nlag"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"lag_separation_in_each_direction", "lag"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"lag_tolerance_in_each_direction", "ltol"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"azimuth", "azm"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"tolerance_of_azimuth", "atol"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"bandwidth_of_azimuth", "bandwh"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"dip_angle", "dip"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"tolerance_of_dip_angle", "dtol"})),
            without_sgems::parse_int_list(get_alias(config, "template_shape", {"dip_bandwidth", "bandwd"})));

        double** vendor_grid = without_sgems::create_vendor_grid_array(grid);
        if(InitialStruct4(key, vendor_grid) == -1)
        {
            throw std::runtime_error("InitialStruct4 failed");
        }

        int** lags = NULL;
        int lag_num = 0;
        if(InitialLagMat(&lags, &lag_num) != 0)
        {
            throw std::runtime_error("InitialLagMat failed");
        }
        if(lag_index < 0 || lag_index >= lag_num)
        {
            throw std::runtime_error("lag index out of range");
        }

        Template* templates = NULL;
        int template_num = 0;
        if(TemplateConstruct_2D(lags[lag_index], pt_Keys->ORDER, &templates, &template_num, 0, pt_Grids->ND) != 0)
        {
            throw std::runtime_error("TemplateConstruct_2D failed");
        }

        const int order = pt_Keys->ORDER;
        double avl = avl_eval(templates->index, order);
        double bell = BellCumulantFromMultiplicityIndex(templates->index, order);
        std::printf("lag_index=%d lag=(%d,%d,%d) template_num=%d root_index=", lag_index, lags[lag_index][0], lags[lag_index][1], lags[lag_index][2], template_num);
        for(int i = 0; i <= order; i++)
        {
            std::printf("%s%d", i == 0 ? "" : ",", templates->index[i]);
        }
        std::printf("\navl=%.17g bell=%.17g diff=%.17g\n", avl, bell, bell - avl);

        const std::vector<std::vector<int>> probes = {
            {1}, {2}, {3}, {4},
            {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4},
            {1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}
        };
        for(const auto& labels : probes)
        {
            std::vector<int> idx;
            if(make_index(order, labels, idx) != 0)
            {
                continue;
            }
            double moment = BellRawMomentFromMultiplicityIndex(idx.data(), order);
            double low = LowOrderCumulant(idx.data(), order);
            double cumulant = BellCumulantFromMultiplicityIndex(idx.data(), order);
            std::printf("labels=");
            for(size_t i = 0; i < labels.size(); i++)
            {
                std::printf("%s%d", i == 0 ? "" : "-", labels[i]);
            }
            std::printf(" index=");
            for(int i = 0; i <= order; i++)
            {
                std::printf("%s%d", i == 0 ? "" : ",", idx[static_cast<size_t>(i)]);
            }
            std::printf(" moment=%.17g low=%.17g bell_cumulant=%.17g low_minus_bell=%.17g\n",
                moment, low, cumulant, low - cumulant);
        }
    }
    catch(const std::exception& e)
    {
        std::fprintf(stderr, "diagnose failed: %s\n", e.what());
        return 1;
    }
    return 0;
}
