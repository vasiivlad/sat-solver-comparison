#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <chrono>

using Clause = std::vector<int>;
using CNF = std::vector<Clause>;

struct Stats {
    double time_sec = 0.0;
    size_t backtracks = 0;
    size_t resolvents = 0;
};

struct Solver {
    CNF formula;
    size_t num_vars = 0;
    std::mt19937 rng;
    Solver() : rng(std::random_device{}()) {}
    void load_dimacs(const std::string& path) {
        std::ifstream in(path);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == 'c') continue;
            if (line[0] == 'p') {
                std::stringstream ss(line);
                std::string tmp;
                ss >> tmp >> tmp >> num_vars;
                continue;
            }
            std::stringstream ss(line);
            Clause cls;
            int lit;
            while (ss >> lit) {
                if (lit == 0) break;
                cls.push_back(lit);
            }
            if (!cls.empty()) {
                std::sort(cls.begin(), cls.end());
                formula.push_back(cls);
            }
        }
    }
    bool any_empty(const CNF& f) const {
        for (auto& c : f) if (c.empty()) return true;
        return false;
    }
    bool solve_resolution(Stats& s) {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<Clause> db = formula;
        bool changed = true;
        while (changed) {
            changed = false;
            size_t n = db.size();
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    for (int lit : db[i]) {
                        if (std::binary_search(db[j].begin(), db[j].end(), -lit)) {
                            Clause res;
                            std::set_union(db[i].begin(), db[i].end(),
                                           db[j].begin(), db[j].end(),
                                           std::back_inserter(res));
                            res.erase(std::remove(res.begin(), res.end(), lit), res.end());
                            res.erase(std::remove(res.begin(), res.end(), -lit), res.end());
                            std::sort(res.begin(), res.end());
                            if (res.empty()) {
                                s.resolvents++;
                                auto end = std::chrono::high_resolution_clock::now();
                                s.time_sec = std::chrono::duration<double>(end - start).count();
                                return false;
                            }
                            if (!std::binary_search(db.begin(), db.end(), res)) {
                                db.push_back(res);
                                changed = true;
                                s.resolvents++;
                            }
                        }
                    }
                }
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        s.time_sec = std::chrono::duration<double>(end - start).count();
        return true;
    }
    bool solve_dp(Stats& s) {
        auto start = std::chrono::high_resolution_clock::now();
        CNF f = formula;
        for (size_t var = 1; var <= num_vars; ++var) {
            CNF pos, neg, rest;
            for (auto& c : f) {
                if (std::binary_search(c.begin(), c.end(), int(var))) pos.push_back(c);
                else if (std::binary_search(c.begin(), c.end(), -int(var))) neg.push_back(c);
                else rest.push_back(c);
            }
            for (auto& a : pos)
                for (auto& b : neg) {
                    Clause res;
                    std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                                   std::back_inserter(res));
                    res.erase(std::remove(res.begin(), res.end(), var), res.end());
                    res.erase(std::remove(res.begin(), res.end(), -int(var)), res.end());
                    std::sort(res.begin(), res.end());
                    if (res.empty()) {
                        s.resolvents++;
                        auto end = std::chrono::high_resolution_clock::now();
                        s.time_sec = std::chrono::duration<double>(end - start).count();
                        return false;
                    }
                    rest.push_back(res);
                    s.resolvents++;
                }
            f.swap(rest);
            if (f.empty()) {
                auto end = std::chrono::high_resolution_clock::now();
                s.time_sec = std::chrono::duration<double>(end - start).count();
                return true;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        s.time_sec = std::chrono::duration<double>(end - start).count();
        return true;
    }
    void unit_propagate(CNF& f, std::vector<int>& asn) {
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto it = f.begin(); it != f.end();) {
                if (it->size() == 1) {
                    int lit = (*it)[0];
                    asn[std::abs(lit)] = lit > 0 ? 1 : -1;
                    for (auto it2 = f.begin(); it2 != f.end();) {
                        if (std::binary_search(it2->begin(), it2->end(), lit)) {
                            it2 = f.erase(it2);
                            continue;
                        }
                        if (std::binary_search(it2->begin(), it2->end(), -lit)) {
                            it2->erase(std::lower_bound(it2->begin(), it2->end(), -lit));
                        }
                        ++it2;
                    }
                    changed = true;
                    break;
                } else ++it;
            }
        }
    }
    int pick_var(const CNF& f, const std::string& h) {
        if (h == "random") {
            std::uniform_int_distribution<int> dist(1, int(num_vars));
            int v;
            do { v = dist(rng); } while (false);
            return v;
        }
        if (h == "moms") {
            size_t kmin = SIZE_MAX;
            for (auto& c : f) kmin = std::min(kmin, c.size());
            std::unordered_map<int, size_t> pos, neg;
            for (auto& c : f) if (c.size() == kmin)
                for (int lit : c) (lit > 0 ? pos[lit] : neg[-lit])++;
            size_t best = 0; int bestv = 1;
            for (size_t v = 1; v <= num_vars; ++v) {
                size_t s = (pos[v] + neg[v]) * ((1u << kmin) - 1) + pos[v] * neg[v];
                if (s > best) { best = s; bestv = int(v); }
            }
            return bestv;
        }
        std::vector<double> jw(num_vars + 1);
        for (auto& c : f) {
            double w = 1.0 / (1u << c.size());
            for (int lit : c)
                jw[std::abs(lit)] += w;
        }
        return int(std::max_element(jw.begin() + 1, jw.end()) - jw.begin());
    }
    bool dpll_rec(CNF& f, std::vector<int>& asn, Stats& s,
                  const std::string& heuristic) {
        if (f.empty()) return true;
        if (any_empty(f)) return false;
        unit_propagate(f, asn);
        if (f.empty()) return true;
        if (any_empty(f)) return false;
        int var = pick_var(f, heuristic);
        if (var == 0) return false;
        for (int sign : {1, -1}) {
            s.backtracks++;
            CNF g = f;
            asn[var] = sign;
            for (auto it = g.begin(); it != g.end();) {
                if (std::binary_search(it->begin(), it->end(), sign * var)) {
                    it = g.erase(it);
                    continue;
                }
                if (std::binary_search(it->begin(), it->end(), -sign * var)) {
                    it->erase(std::lower_bound(it->begin(), it->end(), -sign * var));
                }
                ++it;
            }
            if (dpll_rec(g, asn, s, heuristic)) return true;
            asn[var] = 0;
        }
        return false;
    }
    bool solve_dpll(const std::string& heuristic, Stats& s) {
        auto start = std::chrono::high_resolution_clock::now();
        CNF f = formula;
        std::vector<int> assignment(num_vars + 1, 0);
        bool sat = dpll_rec(f, assignment, s, heuristic);
        auto end = std::chrono::high_resolution_clock::now();
        s.time_sec = std::chrono::duration<double>(end - start).count();
        return sat;
    }
};

#endif
