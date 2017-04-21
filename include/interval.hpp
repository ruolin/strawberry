//
// Created by Ruolin Liu on 11/25/16.
//

#ifndef STRAWBERRY_INTERVAL_RANGES_H
#define STRAWBERRY_INTERVAL_RANGES_H

#include<vector>
#include<cassert>
#include<iostream>

//template <typename TValue, typename TCargo >
//class CargoInterval
//{
//public:
//    TValue i1;
//    TValue i2;
//    TCargo cargo;
//
//    CargoInterval() : i1(), i2(), cargo() {}
//    CargoInterval(TValue i1, TValue i2, TCargo cargo) : i1(i1), i2(i2), cargo(cargo) {}
//};

class BasicInterval{
    /*
     * half open interval encoding.
     * e.g. [left, right) -> coverage
     */
private:
    int left_;
    int right_;
public:
    using TDepth = int ;
    BasicInterval(): left_(0), right_(0){}
    BasicInterval(int l, int r) : left_(l), right_(r) {
        std::cerr<<"interval "<<l<<"-"<<r<<" constructed\n";
    }
    void left(int l) {left_ = l;}
    void right(int r) {right_ = r;}
    decltype(auto) left() const {return left_;}
    decltype(auto) right() const {return right_;}
    int depth() const {return 1;}

    friend std::ostream& operator<<(std::ostream&, const BasicInterval& );
};

inline std::ostream& operator<<(std::ostream& os, const BasicInterval& bi){
    os <<"interval ["<<bi.left() <<"-"<<bi.right()<<")"<<std::endl;
    return os;
}

template<typename TInterval, bool half_open>
inline decltype (auto) coverage(const std::vector<TInterval>& );


template<typename TInterval = BasicInterval , bool half_open = true>
class IRanges{
    /*
     * Use BasicInterval as default interval representation inside class
     */

private:

//    template<typename DataType = int, typename CountType = int>
//    struct RLEncoding {
//        DataType data ;
//        CountType count;
//        //std::make_pair<Type, TLength> duo;
//    };

    std::vector<TInterval> invs_;
//    std::vector<RLEncoding<>> rles_;
    int left_;
    int right_;

    void setBoundries(){
        auto lmin = std::min_element(invs_.begin(), invs_.end(), [&](auto lhs, auto rhs ){
            return lhs.left() < rhs.left();
        });
        auto rmax = std::max_element(invs_.begin(), invs_.end(), [&](auto lhs, auto rhs ){
            return lhs.right() < rhs.right();
        });
        left_ = lmin->left();
        right_ = rmax->right();
    }

public:
    template<typename IntVec = std::initializer_list<int> >
    IRanges(IntVec starts, IntVec ends) {
        assert (starts.size() == ends.size());
        auto e = ends.begin();
        for (auto s = starts.begin(); s < starts.end(); ++s) {
            auto end_pos = *e ;
            if (!half_open) ++end_pos;  // if input is a close interval, convert to half open interval
            TInterval inv(*s, end_pos);
            invs_.push_back(move(inv));
            ++e;
        }
        setBoundries();
    }

    IRanges(std::vector<TInterval> invs): invs_(invs) {
        if (!half_open) {
           for (TInterval& inv: invs_) {
               inv.right(inv.right() + 1);
           }
        }
        setBoundries();
    }

    std::vector<TInterval> reduce(){
        /*
         *  Merge redundant	ranges,	and return the minimum
         *  non-overlapping ranges covering all the input ranges.
         */
        std::vector<typename TInterval::TDepth> cov = coverage<TInterval, true> (invs_);
        std::vector<TInterval> result;

        bool prev_base_is_covered = false;
        bool open = false;

        for (auto it = cov.begin(); it != cov.end(); ++it) {
            if (*it > 0) {
                if (prev_base_is_covered) {

                } else {
                    TInterval inv;
                    inv.left(std::distance(cov.begin(), it));
                    result.push_back(inv);
                    open = true;
                }
                prev_base_is_covered = true;
            } else {
                if (prev_base_is_covered) {
                    result.back().right(std::distance(cov.begin(), it));
                    open = false;
                } else {

                }
                prev_base_is_covered = false;
            }
        }
        if (open) {
            result.back().right(cov.size());
        }

        return result;
    }

    std::vector<TInterval> disjoint(){
    /*
     * return non-overlapping intervals
     */
        std::vector<typename TInterval::TDepth> cov = coverage<TInterval, true> (invs_);
        std::vector<int> bars;
        std::vector<TInterval> result;
        for (auto& inv: invs_) {
            bars.push_back( inv.left() );
            bars.push_back( inv.right() );
        }
        sort(bars.begin(), bars.end());
        auto last = unique(bars.begin(), bars.end());
        bars.erase(last, bars.end());

        bool is_left = true;
        for (auto it = bars.begin(); it != bars.end(); ++it) {
            if (is_left) {
                TInterval sub_inv;
                sub_inv.left(*it);
                result.push_back(sub_inv);
                is_left = false;
            } else {
                if (*it == result.back().left()) {
                    result.pop_back();
                }
                else {
                    result.back().right(*it);
                    if (cov[*it - left_ ] > 0) --it;
                }
                is_left = true;
            }
        }
        if (!is_left) {
            result.pop_back();
        }

        if (!half_open) {
            for (auto & res: result) res.right( res.right() -1);
        }
        return result;
    }
};

template<typename TInterval, bool half_open>
inline decltype (auto) coverage(std::vector<TInterval>& invs){

typedef typename TInterval::TDepth CovType;

auto lmin = std::min_element(invs.begin(), invs.end(), [&](auto lhs, auto rhs ){
   return lhs.left() < rhs.left();
});
auto rmax = std::max_element(invs.begin(), invs.end(), [&](auto lhs, auto rhs ){
    return lhs.right() < rhs.right();
});
int left = lmin->left();
int right = rmax->right();
assert(left < right);

std::vector<CovType> coverage;
if (half_open) coverage.resize(right - left, 0);
else coverage.resize(right - left + 1, 0);
for (auto s = invs.begin(); s != invs.end(); ++s) {
    int end = 0;

    if (half_open) end = s->right();
    else end = s->right() + 1;

        for (int p = s->left(); p < end; ++p) {
             coverage[p-left] += s->depth();
        }
    }
    return coverage;
}

#endif //STRAWBERRY_INTERVAL_RANGES_H
