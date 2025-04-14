#include <bits/stdc++.h>
#pragma GCC optimize("Ofast")

using namespace std;
using namespace chrono;

typedef unsigned long long ull;
typedef pair<ull, ull> p_ull;
typedef pair<int, int> pii;
typedef pair<int, bool> pib;
typedef vector<int> vi;
typedef vector<int>::iterator vi_it;
typedef pair<vi_it, vi_it> p_vi_it;
typedef void (*sort_fcn)(vi_it, vi_it);

const int SEED = 1234;
const int _1K = 1e3;
const int _10K = 1e4;
const int _100K = 1e5;
const int _1M = 1e6;
const int INF = 0x3f3f3f3f;

mt19937 input_data_gen(SEED);
mt19937 quick_sort_gen(SEED);
mt19937 random_shuffle_gen(SEED);

vi temp(_1M*4);

struct Exp_Result {
    int dur[4][4];
    bool is_valid[4][4];
};

struct Experiment {
    const int LOW = -1e9;
    const int HIGH = 1e9;
    const int NLOGN_NUM_REP = 10;
    const int N2_NUM_REP = 2;
    const double PARTIALLY_SORTED_SWAP_RATE = 0.01;

    enum Input_Type {
        SORTED,
        REVERSE_SORTED,
        RAND,
        PARTIALLY_SORTED
    };

    vector<vector<int>> input_data[4];
    uniform_int_distribution<int> dis;
    ofstream out_file;

    Experiment() : dis(LOW, HIGH), out_file("result.txt") {
        gen_sorted();
        gen_reverse_sorted();
        gen_rand();
        gen_partially_sorted();
    }

    p_ull compute_hash(vi_it begin, vi_it end) {
        ull hash_sum = 0;
        ull hash_xor = 0;
        const ull OFFSET = 1e9;

        for(vi_it it = begin; it != end; ++it) {
            ull val = *it + OFFSET;
            hash_sum += val;
            hash_xor ^= val;
        }

        return {hash_sum, hash_xor};
    }

    bool validation_test(vi_it unsorted_begin, vi_it unsorted_end, vi_it sorted_begin, vi_it sorted_end) {
        int sz1 = unsorted_end-unsorted_begin;
        int sz2 = sorted_end-sorted_begin;
        if(sz1 != sz2) return 0;
    
        bool flag = 1;
        for(int i = 1; i < sz1; ++i)
            if(*(sorted_begin+i-1) > *(sorted_begin+i)) {
                flag = 0;
                break;
            }
        if(!flag) return 0;

        p_ull hash_unsorted = compute_hash(unsorted_begin, unsorted_end);
        p_ull hash_sorted = compute_hash(sorted_begin, sorted_end);
        return hash_unsorted == hash_sorted;
    }

    pib nlogn_eval(vi_it unsorted_begin, vi_it unsorted_end, sort_fcn custom_sort) {
        int sz = unsorted_end-unsorted_begin;
        vi sorted(sz);

        int sum = 0;
        for(int i = 0; i < NLOGN_NUM_REP; ++i) {
            vi_it it_sorted = sorted.begin();
            for(vi_it it = unsorted_begin; it != unsorted_end; ++it)
                *it_sorted++ = *it;
            system_clock::time_point start_time = system_clock::now();
            custom_sort(sorted.begin(), sorted.end());
            system_clock::time_point end_time = system_clock::now();
            if(!validation_test(unsorted_begin, unsorted_end, sorted.begin(), sorted.end()))
                return {-1, 0};
            milliseconds ms = duration_cast<milliseconds>(end_time-start_time);
            sum += ms.count();
        }
        return {sum/NLOGN_NUM_REP, 1};
    }

    Exp_Result nlogn_test(sort_fcn custom_sort) {
        Exp_Result res;
        for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j) {
            pib ret = nlogn_eval(input_data[i][j].begin(), input_data[i][j].end(), custom_sort);
            res.dur[i][j] = ret.first;
            res.is_valid[i][j] = ret.second;
        }
        return res;
    }

    pib n2_eval(vi_it unsorted_begin, vi_it unsorted_end, sort_fcn custom_sort) {
        int sz = unsorted_end-unsorted_begin;
        vi sorted(sz);

        int sum = 0;
        for(int i = 0; i < N2_NUM_REP; ++i) {
            vi_it it_sorted = sorted.begin();
            for(vi_it it = unsorted_begin; it != unsorted_end; ++it)
                *it_sorted++ = *it;
            system_clock::time_point start_time = system_clock::now();
            custom_sort(sorted.begin(), sorted.end());
            system_clock::time_point end_time = system_clock::now();
            if(!validation_test(unsorted_begin, unsorted_end, sorted.begin(), sorted.end()))
                return {-1, 0};
            milliseconds ms = duration_cast<milliseconds>(end_time-start_time);
            sum += ms.count();
        }
        return {sum/N2_NUM_REP, 1};
    }

    Exp_Result n2_test(sort_fcn custom_sort) {
        Exp_Result res;
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 3; ++j) {
                pib ret = nlogn_eval(input_data[i][j].begin(), input_data[i][j].end(), custom_sort);
                res.dur[i][j] = ret.first;
                res.is_valid[i][j] = ret.second;
            }
            pib ret = n2_eval(input_data[i][3].begin(), input_data[i][3].end(), custom_sort);
            res.dur[i][3] = ret.first;
            res.is_valid[i][3] = ret.second;
        }
        return res;
    }

    void write_file(Exp_Result res, string name) {
        const int WIDTH = 50;
        string line1, line2;
        string data_type[4] = {"Sorted", "Reverse Sorted", "Random", "Partially Sorted"};
        string data_size[4] = {"1K", "10K", "100K", "1M"};

        for(int i = 0; i < WIDTH; ++i) {
            line1.push_back('-');
            line2.push_back('=');
        }

        int left = (WIDTH-name.length()-2)/2;
        int right = (WIDTH-name.length()-2)/2 + (WIDTH-name.length()-2)%2;
        out_file << line2.substr(0, left)+' ' << name << ' '+line2.substr(0, right) << endl;

        for(int i = 0; i < 4; ++i) {
            out_file << endl;
            out_file << data_type[i] << endl;
            for(int j = 0; j < 4; ++j) {
                out_file << data_size[j] << ": " << res.dur[i][j] << "ms, " << (res.is_valid[i][j] ? "valid" : "invalid") << endl;
            }
            out_file << endl;
            if(i < 3) out_file << line1 << endl;
        }
        
        out_file << line2 << endl;
    }

    void gen_sorted() {
        vector<int> vt1k(_1K), vt10k(_10K), vt100k(_100K), vt1m(_1M);

        for(int i = 0; i < _1K; ++i)
            vt1k[i] = i;
        for(int i = 0; i < _10K; ++i)
            vt10k[i] = i;
        for(int i = 0; i < _100K; ++i)
            vt100k[i] = i;
        for(int i = 0; i < _1M; ++i)
            vt1m[i] = i;
        
        input_data[SORTED].push_back(vt1k);
        input_data[SORTED].push_back(vt10k);
        input_data[SORTED].push_back(vt100k);
        input_data[SORTED].push_back(vt1m);
    }

    void gen_reverse_sorted() {
        vector<int> vt1k(_1K), vt10k(_10K), vt100k(_100K), vt1m(_1M);

        for(int i = 0; i < _1K; ++i)
            vt1k[i] = -i;
        for(int i = 0; i < _10K; ++i)
            vt10k[i] = -i;
        for(int i = 0; i < _100K; ++i)
            vt100k[i] = -i;
        for(int i = 0; i < _1M; ++i)
            vt1m[i] = -i;
        
        input_data[REVERSE_SORTED].push_back(vt1k);
        input_data[REVERSE_SORTED].push_back(vt10k);
        input_data[REVERSE_SORTED].push_back(vt100k);
        input_data[REVERSE_SORTED].push_back(vt1m);
    }

    void gen_rand() {
        vector<int> vt1k(_1K), vt10k(_10K), vt100k(_100K), vt1m(_1M);

        for(int i = 0; i < _1K; ++i)
            vt1k[i] = dis(input_data_gen);
        for(int i = 0; i < _10K; ++i)
            vt10k[i] = dis(input_data_gen);
        for(int i = 0; i < _100K; ++i)
            vt100k[i] = dis(input_data_gen);
        for(int i = 0; i < _1M; ++i)
            vt1m[i] = dis(input_data_gen);
        
        input_data[RAND].push_back(vt1k);
        input_data[RAND].push_back(vt10k);
        input_data[RAND].push_back(vt100k);
        input_data[RAND].push_back(vt1m);
    }

    void gen_partially_sorted() {
        vector<int> vt1k(_1K), vt10k(_10K), vt100k(_100K), vt1m(_1M);
        uniform_int_distribution<int> dis1k(0, _1K);
        uniform_int_distribution<int> dis10k(0, _10K);
        uniform_int_distribution<int> dis100k(0, _100K);
        uniform_int_distribution<int> dis1m(0, _1M);

        for(int i = 0; i < _1K; ++i)
            vt1k[i] = i;
        for(int k = 0; k < _1K*PARTIALLY_SORTED_SWAP_RATE; ++k) {
            int i = dis1k(input_data_gen), j;
            do {
                j = dis1k(input_data_gen);
            } while(i == j);
            swap(vt1k[i], vt1k[j]);
        }

        for(int i = 0; i < _10K; ++i)
            vt10k[i] = i;
        for(int k = 0; k < _10K*PARTIALLY_SORTED_SWAP_RATE; ++k) {
            int i = dis10k(input_data_gen), j;
            do {
                j = dis10k(input_data_gen);
            } while(i == j);
            swap(vt10k[i], vt10k[j]);
        }
        
        for(int i = 0; i < _100K; ++i)
            vt100k[i] = i;
        for(int k = 0; k < _100K*PARTIALLY_SORTED_SWAP_RATE; ++k) {
            int i = dis100k(input_data_gen), j;
            do {
                j = dis100k(input_data_gen);
            } while(i == j);
            swap(vt100k[i], vt100k[j]);
        }

        for(int i = 0; i < _1M; ++i)
            vt1m[i] = i;
        for(int k = 0; k < _1M*PARTIALLY_SORTED_SWAP_RATE; ++k) {
            int i = dis1m(input_data_gen), j;
            do {
                j = dis1m(input_data_gen);
            } while(i == j);
            swap(vt1m[i], vt1m[j]);
        }
        
        input_data[PARTIALLY_SORTED].push_back(vt1k);
        input_data[PARTIALLY_SORTED].push_back(vt10k);
        input_data[PARTIALLY_SORTED].push_back(vt100k);
        input_data[PARTIALLY_SORTED].push_back(vt1m);
    }
};

struct Segment_Tree {
    int tree[1 << 21];
    int sz;

    void init(vi_it begin, vi_it end) {
        int size = end-begin;
        if(size <= _1K) sz = 1 << 10;
        else if(size <= _10K) sz = 1 << 14;
        else if(size <= _100K) sz = 1 << 17;
        else if(size <= _1M) sz = 1 << 20;

        for(int i = sz; i < sz+size; ++i)
            tree[i] = *(begin+i-sz);
        for(int i = sz+size; i < 2*sz; ++i)
            tree[i] = INF;

        for(int i = sz-1; i; --i)
            tree[i] = min(tree[i << 1], tree[i << 1 | 1]);
    }

    int remove_root() {
        int x = 1, ret = tree[1];
        while(x < sz) {
            if(tree[x] == tree[x << 1]) x <<= 1;
            else x = x << 1 | 1;
            tree[x >> 1] = INF;
        }
        tree[x] = INF;

        while(x >>= 1)
            tree[x] = min(tree[x << 1], tree[x << 1 | 1]);

        return ret;
    }
}seg;

struct Frame {
    vi_it begin;
    vi_it end;
    bool sorted;
};

void merge_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    merge_sort(begin, begin+sz/2);
    merge_sort(begin+sz/2, end);

    vi_it it1 = begin, it2 = begin+sz/2;
    vi_it e1 = it2, e2 = end;
    vi_it it_temp = temp.begin();
    
    while(it1 != e1 && it2 != e2)
        if(*it1 <= *it2) *it_temp++ = *it1++;
        else *it_temp++ = *it2++;

    while(it1 != e1) *it_temp++ = *it1++;
    while(it2 != e2) *it_temp++ = *it2++;

    it_temp = temp.begin();
    for(vi_it it = begin; it != end; ++it)
        *it = *it_temp++;
}

void non_rec_merge_sort(vi_it begin, vi_it end) {
    stack<Frame> st;
    st.push({begin, end, 0});

    while(!st.empty()) {
        Frame cur = st.top();
        st.pop();
        vi_it b = cur.begin;
        vi_it e = cur.end;
        bool sorted = cur.sorted;
        int sz = e-b;
        if(sz <= 1) continue;

        if(!sorted) {
            st.push({b, e, 1});
            st.push({b, b+sz/2, 0});
            st.push({b+sz/2, e, 0});
            continue;
        }

        vi_it it1 = b, it2 = b+sz/2;
        vi_it e1 = it2, e2 = e;
        vi_it it_temp = temp.begin();
        
        while(it1 != e1 && it2 != e2)
            if(*it1 <= *it2) *it_temp++ = *it1++;
            else *it_temp++ = *it2++;
            
        while(it1 != e1) *it_temp++ = *it1++;
        while(it2 != e2) *it_temp++ = *it2++;

        it_temp = temp.begin();
        for(vi_it it = b; it != e; ++it)
            *it = *it_temp++;
    }
}

void three_way_merge_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;
    if(sz == 2) {
        merge_sort(begin, end);
        return;
    }

    int gap = sz/3;
    three_way_merge_sort(begin, begin+gap);
    three_way_merge_sort(begin+gap, begin+2*gap);
    three_way_merge_sort(begin+2*gap, end);

    vi_it it1 = begin, it2 = begin+gap, it3 = begin+2*gap;
    vi_it e1 = it2, e2 = it3, e3 = end;
    vi_it it_temp = temp.begin();
    
    while(it1 != e1 && it2 != e2 && it3 != e3)
        if(*it1 <= *it2 && *it1 <= *it3) *it_temp++ = *it1++;
        else if(*it2 <= *it3) *it_temp++ = *it2++;
        else *it_temp++ = *it3++;

    while(it1 != e1 && it2 != e2)
        if(*it1 <= *it2) *it_temp++ = *it1++;
        else *it_temp++ = *it2++;
    while(it2 != e2 && it3 != e3)
        if(*it2 <= *it3) *it_temp++ = *it2++;
        else *it_temp++ = *it3++;
    while(it1 != e1 && it3 != e3)
        if(*it1 <= *it3) *it_temp++ = *it1++;
        else *it_temp++ = *it3++;

    while(it1 != e1) *it_temp++ = *it1++;
    while(it2 != e2) *it_temp++ = *it2++;
    while(it3 != e3) *it_temp++ = *it3++;

    it_temp = temp.begin();
    for(vi_it it = begin; it != end; ++it)
        *it = *it_temp++;
}

void non_rec_three_way_merge_sort(vi_it begin, vi_it end) {
    stack<Frame> st;
    st.push({begin, end, 0});

    while(!st.empty()) {
        Frame cur = st.top();
        st.pop();
        vi_it b = cur.begin;
        vi_it e = cur.end;
        bool sorted = cur.sorted;
        int sz = e-b;
        int gap = sz/3;
        if(sz <= 1) continue;
        if(sz == 2) {
            merge_sort(b, e);
            continue;
        }

        if(!sorted) {
            st.push({b, e, 1});
            st.push({b, b+gap, 0});
            st.push({b+gap, b+2*gap, 0});
            st.push({b+2*gap, e, 0});
            continue;
        }

        vi_it it1 = b, it2 = b+gap, it3 = b+2*gap;
        vi_it e1 = it2, e2 = it3, e3 = e;
        vi_it it_temp = temp.begin();

        while(it1 != e1 && it2 != e2 && it3 != e3)
            if(*it1 <= *it2 && *it1 <= *it3) *it_temp++ = *it1++;
            else if(*it2 <= *it3) *it_temp++ = *it2++;
            else *it_temp++ = *it3++;
        
        while(it1 != e1 && it2 != e2)
            if(*it1 <= *it2) *it_temp++ = *it1++;
            else *it_temp++ = *it2++;
        while(it2 != e2 && it3 != e3)
            if(*it2 <= *it3) *it_temp++ = *it2++;
            else *it_temp++ = *it3++;
        while(it1 != e1 && it3 != e3)
            if(*it1 <= *it3) *it_temp++ = *it1++;
            else *it_temp++ = *it3++;

        while(it1 != e1) *it_temp++ = *it1++;
        while(it2 != e2) *it_temp++ = *it2++;
        while(it3 != e3) *it_temp++ = *it3++;

        it_temp = temp.begin();
        for(vi_it it = b; it != e; ++it)
            *it = *it_temp++;
    }
}

void four_way_merge_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;
    if(sz == 2) {
        merge_sort(begin, end);
        return;
    }
    if(sz == 3) {
        three_way_merge_sort(begin, end);
        return;
    }

    int gap = sz/4;
    four_way_merge_sort(begin, begin+gap);
    four_way_merge_sort(begin+gap, begin+2*gap);
    four_way_merge_sort(begin+2*gap, begin+3*gap);
    four_way_merge_sort(begin+3*gap, end);

    vi_it it1 = begin, it2 = begin+gap, it3 = begin+2*gap, it4 = begin+3*gap;
    vi_it e1 = it2, e2 = it3, e3 = it4, e4 = end;
    vi_it it_temp = temp.begin();
    
    while(it1 != e1 && it2 != e2 && it3 != e3 && it4 != e4)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it4) *it_temp++ = *it2++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;

    while(it1 != e1 && it2 != e2 && it3 != e3)
        if(*it1 <= *it2 && *it1 <= *it3) *it_temp++ = *it1++;
        else if(*it2 <= *it3) *it_temp++ = *it2++;
        else *it_temp++ = *it3++;
    while(it1 != e1 && it2 != e2 && it4 != e4)
        if(*it1 <= *it2 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it2 <= *it4) *it_temp++ = *it2++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it3 != e3 && it4 != e4)
        if(*it1 <= *it3 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it2 != e2 && it3 != e3 && it4 != e4)
        if(*it2 <= *it3 && *it2 <= *it4) *it_temp++ = *it2++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;

    while(it1 != e1 && it2 != e2)
        if(*it1 <= *it2) *it_temp++ = *it1++;
        else *it_temp++ = *it2++;
    while(it1 != e1 && it3 != e3)
        if(*it1 <= *it3) *it_temp++ = *it1++;
        else *it_temp++ = *it3++;
    while(it1 != e1 && it4 != e4)
        if(*it1 <= *it4) *it_temp++ = *it1++;
        else *it_temp++ = *it4++;
    while(it2 != e2 && it3 != e3)
        if(*it2 <= *it3) *it_temp++ = *it2++;
        else *it_temp++ = *it3++;
    while(it2 != e2 && it4 != e4)
        if(*it2 <= *it4) *it_temp++ = *it2++;
        else *it_temp++ = *it4++;
    while(it3 != e3 && it4 != e4)
        if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;

    while(it1 != e1) *it_temp++ = *it1++;
    while(it2 != e2) *it_temp++ = *it2++;
    while(it3 != e3) *it_temp++ = *it3++;
    while(it4 != e4) *it_temp++ = *it4++;

    it_temp = temp.begin();
    for(vi_it it = begin; it != end; ++it)
        *it = *it_temp++;
}

void five_way_merge_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;
    if(sz == 2) {
        merge_sort(begin, end);
        return;
    }
    if(sz == 3) {
        three_way_merge_sort(begin, end);
        return;
    }
    if(sz == 4) {
        four_way_merge_sort(begin, end);
        return;
    }

    int gap = sz/5;
    five_way_merge_sort(begin, begin+gap);
    five_way_merge_sort(begin+gap, begin+2*gap);
    five_way_merge_sort(begin+2*gap, begin+3*gap);
    five_way_merge_sort(begin+3*gap, begin+4*gap);
    five_way_merge_sort(begin+4*gap, end);

    vi_it it1 = begin, it2 = begin+gap, it3 = begin+2*gap, it4 = begin+3*gap, it5 = begin+4*gap;
    vi_it e1 = begin+gap, e2 = begin+2*gap, e3 = begin+3*gap, e4 = begin+4*gap, e5 = end;

    vi_it it_temp = temp.begin();
    while(it1 != e1 && it2 != e2 && it3 != e3 && it4 != e4 && it5 != e5)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it4 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it4 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it3 <= *it4 && *it3 <= *it5) *it_temp++ = *it3++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;

    while(it1 != e1 && it2 != e2 && it3 != e3 && it4 != e4)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it4) *it_temp++ = *it2++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it2 != e2 && it3 != e3 && it5 != e5)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it3 <= *it5) *it_temp++ = *it3++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it2 != e2 && it4 != e4 && it5 != e5)
        if(*it1 <= *it2 && *it1 <= *it4 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it2 <= *it4 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it3 != e3 && it4 != e4 && it5 != e5)
        if(*it1 <= *it3 && *it1 <= *it4 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it3 <= *it4 && *it3 <= *it5) *it_temp++ = *it3++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it2 != e2 && it3 != e3 && it4 != e4 && it5 != e5)
        if(*it2 <= *it3 && *it2 <= *it4 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it3 <= *it4 && *it3 <= *it5) *it_temp++ = *it3++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;

    while(it1 != e1 && it2 != e2 && it3 != e3)
        if(*it1 <= *it2 && *it1 <= *it3) *it_temp++ = *it1++;
        else if(*it2 <= *it3) *it_temp++ = *it2++;
        else *it_temp++ = *it3++;
    while(it1 != e1 && it2 != e2 && it4 != e4)
        if(*it1 <= *it2 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it2 <= *it4) *it_temp++ = *it2++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it2 != e2 && it5 != e5)
        if(*it1 <= *it2 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it2 <= *it5) *it_temp++ = *it2++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it3 != e3 && it4 != e4)
        if(*it1 <= *it3 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it3 != e3 && it5 != e5)
        if(*it1 <= *it3 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it3 <= *it5) *it_temp++ = *it3++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it4 != e4 && it5 != e5)
        if(*it1 <= *it4 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it2 != e2 && it3 != e3 && it4 != e4)
        if(*it2 <= *it3 && *it2 <= *it4) *it_temp++ = *it2++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it2 != e2 && it3 != e3 && it5 != e5)
        if(*it2 <= *it3 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it3 <= *it5) *it_temp++ = *it3++;
        else *it_temp++ = *it5++;
    while(it2 != e2 && it4 != e4 && it5 != e5)
        if(*it2 <= *it4 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it3 != e3 && it4 != e4 && it5 != e5)
        if(*it3 <= *it4 && *it3 <= *it5) *it_temp++ = *it3++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;

    while(it1 != e1 && it2 != e2)
        if(*it1 <= *it2) *it_temp++ = *it1++;
        else *it_temp++ = *it2++;
    while(it1 != e1 && it3 != e3)
        if(*it1 <= *it3) *it_temp++ = *it1++;
        else *it_temp++ = *it3++;
    while(it1 != e1 && it4 != e4)
        if(*it1 <= *it4) *it_temp++ = *it1++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it5 != e5)
        if(*it1 <= *it5) *it_temp++ = *it1++;
        else *it_temp++ = *it5++;
    while(it2 != e2 && it3 != e3)
        if(*it2 <= *it3) *it_temp++ = *it2++;
        else *it_temp++ = *it3++;
    while(it2 != e2 && it4 != e4)
        if(*it2 <= *it4) *it_temp++ = *it2++;
        else *it_temp++ = *it4++;
    while(it2 != e2 && it5 != e5)
        if(*it2 <= *it5) *it_temp++ = *it2++;
        else *it_temp++ = *it5++;
    while(it3 != e3 && it4 != e4)
        if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it3 != e3 && it5 != e5)
        if(*it3 <= *it5) *it_temp++ = *it3++;
        else *it_temp++ = *it5++;
    while(it4 != e4 && it5 != e5)
        if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;

    while(it1 != e1) *it_temp++ = *it1++;
    while(it2 != e2) *it_temp++ = *it2++;
    while(it3 != e3) *it_temp++ = *it3++;
    while(it4 != e4) *it_temp++ = *it4++;
    while(it5 != e5) *it_temp++ = *it5++;

    it_temp = temp.begin();
    for(vi_it it = begin; it != end; ++it)
        *it = *it_temp++;
}

void six_way_merge_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;
    if(sz == 2) {
        merge_sort(begin, end);
        return;
    }
    if(sz == 3) {
        three_way_merge_sort(begin, end);
        return;
    }
    if(sz == 4) {
        four_way_merge_sort(begin, end);
        return;
    }
    if(sz == 5) {
        five_way_merge_sort(begin, end);
        return;
    }

    int gap = sz/6;
    six_way_merge_sort(begin, begin+gap);
    six_way_merge_sort(begin+gap, begin+2*gap);
    six_way_merge_sort(begin+2*gap, begin+3*gap);
    six_way_merge_sort(begin+3*gap, begin+4*gap);
    six_way_merge_sort(begin+4*gap, begin+5*gap);
    six_way_merge_sort(begin+5*gap, end);

    vi_it it1 = begin, it2 = begin+gap, it3 = begin+2*gap, it4 = begin+3*gap, it5 = begin+4*gap, it6 = begin+5*gap;
    vi_it e1 = begin+gap, e2 = begin+2*gap, e3 = begin+3*gap, e4 = begin+4*gap, e5 = begin+5*gap, e6 = end;

    vi_it it_temp = temp.begin();
    while(it1 != e1 && it2 != e2 && it3 != e3 && it4 != e4 && it5 != e5 && it6 != e6)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it4 && *it1 <= *it5 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it4 && *it2 <= *it5 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it3 <= *it4 && *it3 <= *it5 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it4 <= *it5 && *it4 <= *it6) *it_temp++ = *it4++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;

    while(it1 != e1 && it2 != e2 && it3 != e3 && it4 != e4 && it5 != e5)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it4 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it4 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it3 <= *it4 && *it3 <= *it5) *it_temp++ = *it3++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it2 != e2 && it3 != e3 && it4 != e4 && it6 != e6)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it4 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it4 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it3 <= *it4 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it4 <= *it6) *it_temp++ = *it4++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it2 != e2 && it3 != e3 && it5 != e5 && it6 != e6)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it5 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it5 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it3 <= *it5 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it2 != e2 && it4 != e4 && it5 != e5 && it6 != e6)
        if(*it1 <= *it2 && *it1 <= *it4 && *it1 <= *it5 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it2 <= *it4 && *it2 <= *it5 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it4 <= *it5 && *it4 <= *it6) *it_temp++ = *it4++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it3 != e3 && it4 != e4 && it5 != e5 && it6 != e6)
        if(*it1 <= *it3 && *it1 <= *it4 && *it1 <= *it5 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it3 <= *it4 && *it3 <= *it5 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it4 <= *it5 && *it4 <= *it6) *it_temp++ = *it4++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it2 != e2 && it3 != e3 && it4 != e4 && it5 != e5 && it6 != e6)
        if(*it2 <= *it3 && *it2 <= *it4 && *it2 <= *it5 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it3 <= *it4 && *it3 <= *it5 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it4 <= *it5 && *it4 <= *it6) *it_temp++ = *it4++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;

    while(it1 != e1 && it2 != e2 && it3 != e3 && it4 != e4)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it4) *it_temp++ = *it2++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it2 != e2 && it3 != e3 && it5 != e5)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it3 <= *it5) *it_temp++ = *it3++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it2 != e2 && it3 != e3 && it6 != e6)
        if(*it1 <= *it2 && *it1 <= *it3 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it2 <= *it3 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it3 <= *it6) *it_temp++ = *it3++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it2 != e2 && it4 != e4 && it5 != e5)
        if(*it1 <= *it2 && *it1 <= *it4 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it2 <= *it4 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it2 != e2 && it4 != e4 && it6 != e6)
        if(*it1 <= *it2 && *it1 <= *it4 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it2 <= *it4 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it4 <= *it6) *it_temp++ = *it4++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it2 != e2 && it5 != e5 && it6 != e6)
        if(*it1 <= *it2 && *it1 <= *it5 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it2 <= *it5 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it3 != e3 && it4 != e4 && it5 != e5)
        if(*it1 <= *it3 && *it1 <= *it4 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it3 <= *it4 && *it3 <= *it5) *it_temp++ = *it3++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it3 != e3 && it4 != e4 && it6 != e6)
        if(*it1 <= *it3 && *it1 <= *it4 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it3 <= *it4 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it4 <= *it6) *it_temp++ = *it4++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it3 != e3 && it5 != e5 && it6 != e6)
        if(*it1 <= *it3 && *it1 <= *it5 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it3 <= *it5 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it4 != e4 && it5 != e5 && it6 != e6)
        if(*it1 <= *it4 && *it1 <= *it5 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it4 <= *it5 && *it4 <= *it6) *it_temp++ = *it4++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it2 != e2 && it3 != e3 && it4 != e4 && it5 != e5)
        if(*it2 <= *it3 && *it2 <= *it4 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it3 <= *it4 && *it3 <= *it5) *it_temp++ = *it3++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it2 != e2 && it3 != e3 && it4 != e4 && it6 != e6)
        if(*it2 <= *it3 && *it2 <= *it4 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it3 <= *it4 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it4 <= *it6) *it_temp++ = *it4++;
        else *it_temp++ = *it6++;
    while(it2 != e2 && it3 != e3 && it5 != e5 && it6 != e6)
        if(*it2 <= *it3 && *it2 <= *it5 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it3 <= *it5 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it2 != e2 && it4 != e4 && it5 != e5 && it6 != e6)
        if(*it2 <= *it4 && *it2 <= *it5 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it4 <= *it5 && *it4 <= *it6) *it_temp++ = *it4++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it3 != e3 && it4 != e4 && it5 != e5 && it6 != e6)
        if(*it3 <= *it4 && *it3 <= *it5 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it4 <= *it5 && *it4 <= *it6) *it_temp++ = *it4++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;

    while(it1 != e1 && it2 != e2 && it3 != e3)
        if(*it1 <= *it2 && *it1 <= *it3) *it_temp++ = *it1++;
        else if(*it2 <= *it3) *it_temp++ = *it2++;
        else *it_temp++ = *it3++;
    while(it1 != e1 && it2 != e2 && it4 != e4)
        if(*it1 <= *it2 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it2 <= *it4) *it_temp++ = *it2++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it2 != e2 && it5 != e5)
        if(*it1 <= *it2 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it2 <= *it5) *it_temp++ = *it2++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it2 != e2 && it6 != e6)
        if(*it1 <= *it2 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it2 <= *it6) *it_temp++ = *it2++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it3 != e3 && it4 != e4)
        if(*it1 <= *it3 && *it1 <= *it4) *it_temp++ = *it1++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it3 != e3 && it5 != e5)
        if(*it1 <= *it3 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it3 <= *it5) *it_temp++ = *it3++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it3 != e3 && it6 != e6)
        if(*it1 <= *it3 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it3 <= *it6) *it_temp++ = *it3++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it4 != e4 && it5 != e5)
        if(*it1 <= *it4 && *it1 <= *it5) *it_temp++ = *it1++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it4 != e4 && it6 != e6)
        if(*it1 <= *it4 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it4 <= *it6) *it_temp++ = *it4++;
        else *it_temp++ = *it6++;
    while(it1 != e1 && it5 != e5 && it6 != e6)
        if(*it1 <= *it5 && *it1 <= *it6) *it_temp++ = *it1++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it2 != e2 && it3 != e3 && it4 != e4)
        if(*it2 <= *it3 && *it2 <= *it4) *it_temp++ = *it2++;
        else if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it2 != e2 && it3 != e3 && it5 != e5)
        if(*it2 <= *it3 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it3 <= *it5) *it_temp++ = *it3++;
        else *it_temp++ = *it5++;
    while(it2 != e2 && it3 != e3 && it6 != e6)
        if(*it2 <= *it3 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it3 <= *it6) *it_temp++ = *it3++;
        else *it_temp++ = *it6++;
    while(it2 != e2 && it4 != e4 && it5 != e5)
        if(*it2 <= *it4 && *it2 <= *it5) *it_temp++ = *it2++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it2 != e2 && it4 != e4 && it6 != e6)
        if(*it2 <= *it4 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it4 <= *it6) *it_temp++ = *it4++;
        else *it_temp++ = *it6++;
    while(it2 != e2 && it5 != e5 && it6 != e6)
        if(*it2 <= *it5 && *it2 <= *it6) *it_temp++ = *it2++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it3 != e3 && it4 != e4 && it5 != e5)
        if(*it3 <= *it4 && *it3 <= *it5) *it_temp++ = *it3++;
        else if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it3 != e3 && it4 != e4 && it6 != e6)
        if(*it3 <= *it4 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it4 <= *it6) *it_temp++ = *it4++;
        else *it_temp++ = *it6++;
    while(it3 != e3 && it5 != e5 && it6 != e6)
        if(*it3 <= *it5 && *it3 <= *it6) *it_temp++ = *it3++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    while(it4 != e4 && it5 != e5 && it6 != e6)
        if(*it4 <= *it5 && *it4 <= *it6) *it_temp++ = *it4++;
        else if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    
    while(it1 != e1 && it2 != e2)
        if(*it1 <= *it2) *it_temp++ = *it1++;
        else *it_temp++ = *it2++;
    while(it1 != e1 && it3 != e3)
        if(*it1 <= *it3) *it_temp++ = *it1++;
        else *it_temp++ = *it3++;
    while(it1 != e1 && it4 != e4)
        if(*it1 <= *it4) *it_temp++ = *it1++;
        else *it_temp++ = *it4++;
    while(it1 != e1 && it5 != e5)
        if(*it1 <= *it5) *it_temp++ = *it1++;
        else *it_temp++ = *it5++;
    while(it1 != e1 && it6 != e6)
        if(*it1 <= *it6) *it_temp++ = *it1++;
        else *it_temp++ = *it6++;
    while(it2 != e2 && it3 != e3)
        if(*it2 <= *it3) *it_temp++ = *it2++;
        else *it_temp++ = *it3++;
    while(it2 != e2 && it4 != e4)
        if(*it2 <= *it4) *it_temp++ = *it2++;
        else *it_temp++ = *it4++;
    while(it2 != e2 && it5 != e5)
        if(*it2 <= *it5) *it_temp++ = *it2++;
        else *it_temp++ = *it5++;
    while(it2 != e2 && it6 != e6)
        if(*it2 <= *it6) *it_temp++ = *it2++;
        else *it_temp++ = *it6++;
    while(it3 != e3 && it4 != e4)
        if(*it3 <= *it4) *it_temp++ = *it3++;
        else *it_temp++ = *it4++;
    while(it3 != e3 && it5 != e5)
        if(*it3 <= *it5) *it_temp++ = *it3++;
        else *it_temp++ = *it5++;
    while(it3 != e3 && it6 != e6)
        if(*it3 <= *it6) *it_temp++ = *it3++;
        else *it_temp++ = *it6++;
    while(it4 != e4 && it5 != e5)
        if(*it4 <= *it5) *it_temp++ = *it4++;
        else *it_temp++ = *it5++;
    while(it4 != e4 && it6 != e6)
        if(*it4 <= *it6) *it_temp++ = *it4++;
        else *it_temp++ = *it6++;
    while(it5 != e5 && it6 != e6)
        if(*it5 <= *it6) *it_temp++ = *it5++;
        else *it_temp++ = *it6++;
    
    while(it1 != e1) *it_temp++ = *it1++;
    while(it2 != e2) *it_temp++ = *it2++;
    while(it3 != e3) *it_temp++ = *it3++;
    while(it4 != e4) *it_temp++ = *it4++;
    while(it5 != e5) *it_temp++ = *it5++;
    while(it6 != e6) *it_temp++ = *it6++;
    
    it_temp = temp.begin();
    for(vi_it it = begin; it != end; ++it)
        *it = *it_temp++;
}

void base0_max_heapify(vi_it begin, vi_it end, int idx) {
    int sz = end-begin;
    if(sz <= 1) return;

    int left = 2*idx+1;
    int right = 2*idx+2;
    int largest = idx;

    if(left < sz && *(begin+left) > *(begin+largest))
        largest = left;
    if(right < sz && *(begin+right) > *(begin+largest))
        largest = right;
    
    if(largest != idx) {
        swap(*(begin+idx), *(begin+largest));
        base0_max_heapify(begin, end, largest);
    }
}

void base0_heap_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    int idx = (sz-1)/2-(sz%2);
    for(; idx >= 0; --idx)
        base0_max_heapify(begin, end, idx);

    for(vi_it it = end-1; it != begin; --it) {
        swap(*begin, *it);
        base0_max_heapify(begin, it, 0);
    }
}

void base0_non_rec_max_heapify(vi_it begin, vi_it end, int idx) {
    int sz = end-begin;
    if(sz <= 1) return;

    int left, right, largest = idx;
    while(1) {
        left = 2*idx+1;
        right = 2*idx+2;

        if(left < sz && *(begin+left) > *(begin+largest))
            largest = left;
        if(right < sz && *(begin+right) > *(begin+largest))
            largest = right;
        
        if(largest != idx) {
            swap(*(begin+idx), *(begin+largest));
            idx = largest;
        }
        else break;
    }
}

void base0_non_rec_heap_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    int idx = (sz-1)/2-(sz%2);
    for(; idx >= 0; --idx)
        base0_non_rec_max_heapify(begin, end, idx);

    for(vi_it it = end-1; it != begin; --it) {
        swap(*begin, *it);
        base0_non_rec_max_heapify(begin, it, 0);
    }
}

void base1_max_heapify(vi_it begin, vi_it end, int idx) {
    int sz = end-begin;
    if(sz <= 1) return;

    int left = idx << 1;
    int right = idx << 1 | 1;
    int largest = idx;

    if(left <= sz && *(begin+left-1) > *(begin+largest-1))
        largest = left;
    if(right <= sz && *(begin+right-1) > *(begin+largest-1))
        largest = right;
    
    if(largest != idx) {
        swap(*(begin+idx-1), *(begin+largest-1));
        base1_max_heapify(begin, end, largest);
    }
}

void base1_heap_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    int idx = sz >> 1;
    for(; idx; --idx)
        base1_max_heapify(begin, end, idx);

    for(vi_it it = end-1; it != begin; --it) {
        swap(*begin, *it);
        base1_max_heapify(begin, it, 1);
    }
}

void base1_non_rec_max_heapify(vi_it begin, vi_it end, int idx) {
    int sz = end-begin;
    if(sz <= 1) return;

    int left, right, largest = idx;
    while(1) {
        left = idx << 1;
        right = idx << 1 | 1;

        if(left <= sz && *(begin+left-1) > *(begin+largest-1))
            largest = left;
        if(right <= sz && *(begin+right-1) > *(begin+largest-1))
            largest = right;
        
        if(largest != idx) {
            swap(*(begin+idx-1), *(begin+largest-1));
            idx = largest;
        }
        else break;
    }
}

void base1_non_rec_heap_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    int idx = sz >> 1;
    for(; idx; --idx)
        base1_non_rec_max_heapify(begin, end, idx);

    for(vi_it it = end-1; it != begin; --it) {
        swap(*begin, *it);
        base1_non_rec_max_heapify(begin, it, 1);
    }
}

void bubble_sort(vi_it begin, vi_it end) {
    while(begin+1 < end) {
        end--;
        for(vi_it it = begin; it != end; ++it)
            if(*it > *(it+1)) swap(*it, *(it+1));
    }
}

void insertion_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    for(vi_it it = begin+1; it != end; ++it) {
        vi_it temp_it = it-1;
        do {
            if(*temp_it > *(temp_it+1)) swap(*temp_it, *(temp_it+1));
            else break;
        } while(temp_it-- != begin);
    }
}

void binary_insertion_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    vi_it e = end;
    for(vi_it it = begin+1; it != end; ++it)
        if(*(it-1) <= *it) {
            e = it;
            break;
        }
    reverse(begin, e);

    for(vi_it it = e; it != end; ++it) {
        int val = *it;
        if(val >= *(it-1)) continue;
        vi_it ins_it = upper_bound(begin, it, val);
        for(vi_it temp_it = it; temp_it != ins_it; --temp_it) *temp_it = *(temp_it-1);
        *ins_it = val;
    }
}

void selection_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;
    vi temp(sz);

    vi_it it_temp = temp.begin();

    for(vi_it it = begin; it != end; ++it)
        *it_temp++ = *it;

    for(vi_it ins_it = begin; ins_it != end; ++ins_it) {
        vi_it mn = temp.begin();
        vi_it temp_end = temp.end();
        for(vi_it it = mn+1; it != temp_end; ++it)
            mn = *mn >= *it ? it : mn;
        *ins_it = *mn;
        temp.erase(mn);
    }
}

vi_it partition(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz == 1) return begin;

    uniform_int_distribution<int> dis(0, sz-1);
    swap(*(end-1), *(begin+dis(quick_sort_gen)));

    int pivot = *(end-1);
    vi_it i = begin;
    for(vi_it j = begin; j != end-1; ++j)
        if(*j <= pivot)
            swap(*i++, *j);
    swap(*i, *(end-1));
    return i;
}

void quick_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    vi_it mid = partition(begin, end);
    quick_sort(begin, mid);
    quick_sort(mid+1, end);
}

void library_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    int range = 1;
    vector<pii> available;
    shuffle(begin, end, random_shuffle_gen);

    while(range < sz) {
        available.clear();
        range <<= 1;

        for(int i = 0; i < range; ++i) {
            if(i & 1 || (i >> 1) >= sz) temp[i] = INF;
            else {
                temp[i] = *(begin+(i >> 1));
                available.push_back({temp[i], i});
            }
        }

        int end_idx = min(range, sz);
        for(int i = range >> 1; i < end_idx; ++i) {
            int val = *(begin+i);
            int s = 0, e = available.size()-1, ins_pos = range-1, ava_idx = available.size();
            while(s <= e) {
                int mid = (s+e) >> 1;
                if(available[mid].first >= val) {
                    ins_pos = available[mid].second;
                    ava_idx = mid;
                    e = mid-1;
                }
                else s = mid+1;
            }

            int offset = 0;
            while(1) {
                if(ins_pos+offset >= range) {
                    offset = -offset;
                    while(temp[ins_pos+offset] != INF) offset--;
                    break;
                }
                if(temp[ins_pos+offset] == INF) break;
                if(ins_pos-offset < 0) {
                    while(temp[ins_pos+offset] != INF) offset++;
                    break;
                }
                if(temp[ins_pos-offset] == INF) {
                    offset = -offset;
                    break;
                }
                offset++;
            }

            int asz = available.size();
            int ins_idx = ava_idx;
            temp[ins_pos+offset] = val;
            if(offset > 0) {
                for(int j = ins_pos+offset; j > ins_pos; --j) {
                    swap(temp[j], temp[j-1]);
                    available[ava_idx].second++;
                    ava_idx++;
                }
                available.insert(available.begin()+ins_idx, {val, ins_pos});
            }
            else if(offset < 0) {
                if(ava_idx == asz) {
                    for(int j = ins_pos+offset; j < ins_pos; ++j) {
                        swap(temp[j], temp[j+1]);
                        ava_idx--;
                        available[ava_idx].second--;
                    }
                    available.insert(available.begin()+ins_idx, {val, range-1});
                }
                else {
                    for(int j = ins_pos+offset; j < ins_pos-1; ++j) {
                        swap(temp[j], temp[j+1]);
                        ava_idx--;
                        available[ava_idx].second--;
                    }
                    available.insert(available.begin()+ins_idx, {val, ins_pos-1});
                }
            }
            else available.push_back({val, range-1});
        }

        vi_it it = begin;
        for(int i = 0; i < range; ++i) {
            if(temp[i] == INF) continue;
            *it++ = temp[i];
        }
    }
}

int get_minrun(int sz) {
    int rem = 0;
    while(sz >= 32) {
        rem |= sz & 1;
        sz >>= 1;
    }
    return sz+rem;
}

void rev_binary_insertion_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    vi_it e = end;
    for(vi_it it = begin+1; it != end; ++it)
        if(*(it-1) < *it) {
            e = it;
            break;
        }
    reverse(begin, e);

    for(vi_it it = e; it != end; ++it) {
        int val = *it;
        if(val > *(it-1)) continue;
        vi_it ins_it = lower_bound(begin, it, val);
        for(vi_it temp_it = it; temp_it != ins_it; --temp_it) *temp_it = *(temp_it-1);
        *ins_it = val;
    }
}

void tim_sort_merge(p_vi_it itv1, p_vi_it itv2) {
    vi_it begin = upper_bound(itv1.first, itv1.second, *itv2.first);
    vi_it end = lower_bound(itv2.first, itv2.second, *(itv1.second-1));
    vi_it it_temp = temp.begin();
    
    int sz1 = itv1.second-begin;
    int sz2 = end-itv2.first;
    int cnt1 = 0, cnt2 = 0;

    if(sz1 <= sz2) {
        for(vi_it it = begin; it != itv1.second; ++it)
            *it_temp++ = *it;

        int idx1 = 0, idx2 = 0;
        vi_it it = begin;
        while(idx1 < sz1 && idx2 < sz2) {
            if(temp[idx1] <= *(itv2.first+idx2)) {
                cnt1++;
                cnt2 = 0;
                if(cnt1 == 3) {
                    int gap = 2, val = *(itv2.first+idx2), prev_idx = idx1;
                    while(idx1+gap < sz1 && temp[idx1+gap] <= val) {
                        idx1 += gap;
                        gap <<= 1;
                    }

                    int s = idx1, e = min(idx1+gap, sz1)-1;
                    while(s <= e) {
                        int mid = (s + e) >> 1;
                        if(temp[mid] <= val) {
                            idx1 = mid;
                            s = mid+1;
                        }
                        else e = mid-1;
                    }
                    idx1++;

                    for(int i = prev_idx; i < idx1; ++i)
                        *it++ = temp[i];
                    cnt1 = 0;
                }
                else *it++ = temp[idx1++];
            }
            else {
                cnt2++;
                cnt1 = 0;
                if(cnt2 == 3) {
                    int gap = 2, val = temp[idx1], prev_idx = idx2;
                    while(idx2+gap < sz2 && *(itv2.first+idx2+gap) < val) {
                        idx2 += gap;
                        gap <<= 1;
                    }
                    
                    int s = idx2, e = min(idx2+gap, sz2)-1;
                    while(s <= e) {
                        int mid = (s + e) >> 1;
                        if(*(itv2.first+mid) < val) {
                            idx2 = mid;
                            s = mid+1;
                        }
                        else e = mid-1;
                    }
                    idx2++;

                    for(int i = prev_idx; i < idx2; ++i)
                        *it++ = *(itv2.first+i);
                    cnt2 = 0;
                }
                else *it++ = *(itv2.first+idx2++);
            }
        }

        while(idx1 < sz1) *it++ = temp[idx1++];
        while(idx2 < sz2) *it++ = *(itv2.first+idx2++);
    }

    else {
        for(vi_it it = itv2.first; it != end; ++it)
            *it_temp++ = *it;

        int idx1 = sz1-1, idx2 = sz2-1;
        vi_it it = end-1;
        while(idx1 >= 0 && idx2 >= 0) {
            if(*(begin+idx1) <= temp[idx2]) {
                cnt2++;
                cnt1 = 0;
                if(cnt2 == 3) {
                    int gap = 2, val = *(begin+idx1), prev_idx = idx2;
                    while(idx2-gap >= 0 && temp[idx2-gap] >= val) {
                        idx2 -= gap;
                        gap <<= 1;
                    }
                    
                    int s = max(idx2-gap, 0), e = idx2;
                    while(s <= e) {
                        int mid = (s + e) >> 1;
                        if(temp[mid] >= val) {
                            idx2 = mid;
                            e = mid-1;
                        }
                        else s = mid+1;
                    }
                    idx2--;

                    for(int i = prev_idx; i > idx2; --i)
                        *it-- = temp[i];
                    cnt2 = 0;
                }
                else *it-- = temp[idx2--];
            }
            else {
                cnt1++;
                cnt2 = 0;
                if(cnt1 == 3) {
                    int gap = 2, val = temp[idx2], prev_idx = idx1;
                    while(idx1-gap >= 0 && *(begin+idx1-gap) > val) {
                        idx1 -= gap;
                        gap <<= 1;
                    }

                    int s = max(idx1-gap, 0), e = idx1;
                    while(s <= e) {
                        int mid = (s + e) >> 1;
                        if(*(begin+mid) > val) {
                            idx1 = mid;
                            e = mid-1;
                        }
                        else s = mid+1;
                    }
                    idx1--;

                    for(int i = prev_idx; i > idx1; --i)
                        *it-- = *(begin+i);
                    cnt1 = 0;
                }
                else *it-- = *(begin+idx1--);
            }
        }

        while(idx1 >= 0) *it-- = *(begin+idx1--);
        while(idx2 >= 0) *it-- = temp[idx2--];
    }
}

void tim_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;
    if(sz <= 32) {
        binary_insertion_sort(begin, end);
        return;
    }

    stack<p_vi_it> st;
    int minrun = get_minrun(sz);
    vi_it b = begin, e;
    do {
        int run = min(minrun, int(end-b));
        e = b+run;
        if(run != 1) {
            if(*b <= *(b+1)) {
                binary_insertion_sort(b, e);
                while(e != end && *(e-1) <= *e) e++;
            }
            else {
                reverse(b, e);
                rev_binary_insertion_sort(b, e);
                reverse(b, e);
                while(e != end && *(e-1) > *e) e++;
                reverse(b, e);
            }
        }

        st.push({b, e});
        bool flag = 0;

        while(st.size() >= 4) {
            p_vi_it top1 = st.top();
            int sz1 = top1.second-top1.first;
            st.pop();
            p_vi_it top2 = st.top();
            int sz2 = top2.second-top2.first;
            st.pop();
            p_vi_it top3 = st.top();
            int sz3 = top3.second-top3.first;
            st.pop();
            p_vi_it top4 = st.top();
            int sz4 = top4.second-top4.first;
            st.pop();
            if(sz4 > sz3+sz2 && sz3 > sz2+sz1 && sz2 > sz1) {
                st.push(top4);
                st.push(top3);
                st.push(top2);
                st.push(top1);
                flag = 1;
                break;
            }
            if(sz4 <= sz3+sz2 || sz3 <= sz2) {
                if(sz4 > sz2) {
                    tim_sort_merge(top3, top2);
                    st.push(top4);
                    st.push({top3.first, top2.second});
                }
                else {
                    tim_sort_merge(top4, top3);
                    st.push({top4.first, top3.second});
                    st.push(top2);
                }
                st.push(top1);
            }
            else {
                st.push(top4);
                if(sz3 > sz1) {
                    tim_sort_merge(top2, top1);
                    st.push(top3);
                    st.push({top2.first, top1.second});
                }
                else {
                    tim_sort_merge(top3, top2);
                    st.push({top3.first, top2.second});
                    st.push(top1);
                }
            }
        }

        while(!flag && st.size() >= 3) {
            p_vi_it top1 = st.top();
            int sz1 = top1.second-top1.first;
            st.pop();
            p_vi_it top2 = st.top();
            int sz2 = top2.second-top2.first;
            st.pop();
            p_vi_it top3 = st.top();
            int sz3 = top3.second-top3.first;
            st.pop();
            if(sz3 > sz2+sz1 && sz2 > sz1) {
                st.push(top3);
                st.push(top2);
                st.push(top1);
                flag = 1;
                break;
            }
            if(sz3 > sz1) {
                tim_sort_merge(top2, top1);
                st.push(top3);
                st.push({top2.first, top1.second});
            }
            else {
                tim_sort_merge(top3, top2);
                st.push({top3.first, top2.second});
                st.push(top1);
            }
        }

        while(!flag && st.size() >= 2) {
            p_vi_it top1 = st.top();
            int sz1 = top1.second-top1.first;
            st.pop();
            p_vi_it top2 = st.top();
            int sz2 = top2.second-top2.first;
            st.pop();
            if(sz2 > sz1) {
                st.push(top2);
                st.push(top1);
                break;
            }
            else {
                tim_sort_merge(top2, top1);
                st.push({top2.first, top1.second});
            }
        }

        b = e;
    } while(e != end);

    while(st.size() > 1) {
        p_vi_it top1 = st.top();
        st.pop();
        p_vi_it top2 = st.top();
        st.pop();
        tim_sort_merge(top2, top1);
        st.push({top2.first, top1.second});
    }
}

void no_galloping_tim_sort_merge(p_vi_it itv1, p_vi_it itv2) {
    vi_it begin = upper_bound(itv1.first, itv1.second, *itv2.first);
    vi_it end = lower_bound(itv2.first, itv2.second, *(itv1.second-1));
    vi_it it_temp = temp.begin();
    
    int sz1 = itv1.second-begin;
    int sz2 = end-itv2.first;

    if(sz1 <= sz2) {
        for(vi_it it = begin; it != itv1.second; ++it)
            *it_temp++ = *it;

        int idx1 = 0, idx2 = 0;
        vi_it it = begin;
        while(idx1 < sz1 && idx2 < sz2) {
            if(temp[idx1] <= *(itv2.first+idx2)) *it++ = temp[idx1++];
            else *it++ = *(itv2.first+idx2++);
        }

        while(idx1 < sz1) *it++ = temp[idx1++];
        while(idx2 < sz2) *it++ = *(itv2.first+idx2++);
    }

    else {
        for(vi_it it = itv2.first; it != end; ++it)
            *it_temp++ = *it;

        int idx1 = sz1-1, idx2 = sz2-1;
        vi_it it = end-1;
        while(idx1 >= 0 && idx2 >= 0) {
            if(*(begin+idx1) <= temp[idx2]) *it-- = temp[idx2--];
            else *it-- = *(begin+idx1--);
        }

        while(idx1 >= 0) *it-- = *(begin+idx1--);
        while(idx2 >= 0) *it-- = temp[idx2--];
    }
}

void no_galloping_tim_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;
    if(sz <= 32) {
        binary_insertion_sort(begin, end);
        return;
    }

    stack<p_vi_it> st;
    int minrun = get_minrun(sz);
    vi_it b = begin, e;
    do {
        int run = min(minrun, int(end-b));
        e = b+run;
        if(run != 1) {
            if(*b <= *(b+1)) {
                binary_insertion_sort(b, e);
                while(e != end && *(e-1) <= *e) e++;
            }
            else {
                reverse(b, e);
                rev_binary_insertion_sort(b, e);
                reverse(b, e);
                while(e != end && *(e-1) > *e) e++;
                reverse(b, e);
            }
        }

        st.push({b, e});
        bool flag = 0;

        while(st.size() >= 4) {
            p_vi_it top1 = st.top();
            int sz1 = top1.second-top1.first;
            st.pop();
            p_vi_it top2 = st.top();
            int sz2 = top2.second-top2.first;
            st.pop();
            p_vi_it top3 = st.top();
            int sz3 = top3.second-top3.first;
            st.pop();
            p_vi_it top4 = st.top();
            int sz4 = top4.second-top4.first;
            st.pop();
            if(sz4 > sz3+sz2 && sz3 > sz2+sz1 && sz2 > sz1) {
                st.push(top4);
                st.push(top3);
                st.push(top2);
                st.push(top1);
                flag = 1;
                break;
            }
            if(sz4 <= sz3+sz2 || sz3 <= sz2) {
                if(sz4 > sz2) {
                    no_galloping_tim_sort_merge(top3, top2);
                    st.push(top4);
                    st.push({top3.first, top2.second});
                }
                else {
                    no_galloping_tim_sort_merge(top4, top3);
                    st.push({top4.first, top3.second});
                    st.push(top2);
                }
                st.push(top1);
            }
            else {
                st.push(top4);
                if(sz3 > sz1) {
                    no_galloping_tim_sort_merge(top2, top1);
                    st.push(top3);
                    st.push({top2.first, top1.second});
                }
                else {
                    no_galloping_tim_sort_merge(top3, top2);
                    st.push({top3.first, top2.second});
                    st.push(top1);
                }
            }
        }

        while(!flag && st.size() >= 3) {
            p_vi_it top1 = st.top();
            int sz1 = top1.second-top1.first;
            st.pop();
            p_vi_it top2 = st.top();
            int sz2 = top2.second-top2.first;
            st.pop();
            p_vi_it top3 = st.top();
            int sz3 = top3.second-top3.first;
            st.pop();
            if(sz3 > sz2+sz1 && sz2 > sz1) {
                st.push(top3);
                st.push(top2);
                st.push(top1);
                flag = 1;
                break;
            }
            if(sz3 > sz1) {
                no_galloping_tim_sort_merge(top2, top1);
                st.push(top3);
                st.push({top2.first, top1.second});
            }
            else {
                no_galloping_tim_sort_merge(top3, top2);
                st.push({top3.first, top2.second});
                st.push(top1);
            }
        }

        while(!flag && st.size() >= 2) {
            p_vi_it top1 = st.top();
            int sz1 = top1.second-top1.first;
            st.pop();
            p_vi_it top2 = st.top();
            int sz2 = top2.second-top2.first;
            st.pop();
            if(sz2 > sz1) {
                st.push(top2);
                st.push(top1);
                break;
            }
            else {
                no_galloping_tim_sort_merge(top2, top1);
                st.push({top2.first, top1.second});
            }
        }

        b = e;
    } while(e != end);

    while(st.size() > 1) {
        p_vi_it top1 = st.top();
        st.pop();
        p_vi_it top2 = st.top();
        st.pop();
        no_galloping_tim_sort_merge(top2, top1);
        st.push({top2.first, top1.second});
    }
}

void cocktail_shaker_sort(vi_it begin, vi_it end) {
    bool flag = 1;
    while(begin+1 < end) {
        if(flag) {
            end--;
            for(vi_it it = begin; it != end; ++it)
                if(*it > *(it+1)) swap(*it, *(it+1));
        }
        else {
            for(vi_it it = end-1; it != begin; --it)
                if(*(it-1) > *it) swap(*(it-1), *it);
            begin++;
        }
        flag = !flag;
    }
}

void comb_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;
    
    int gap = end-begin;
    double shrink = 1.3;
    bool sorted = 0;

    while(!sorted) {
        gap = floor(gap/shrink);
        if(gap <= 1) {
            gap = 1;
            sorted = 1;
        }

        for(int i = 0; i+gap < sz; ++i) {
            if(*(begin+i) > *(begin+i+gap)) {
                swap(*(begin+i), *(begin+i+gap));
                sorted = 0;
            }
        }
    }
}

void tournament_sort(vi_it begin, vi_it end) {
    seg.init(begin, end);
    
    for(vi_it it = begin; it != end; ++it)
        *it = seg.remove_root();
}

void intro_sort_rec(vi_it begin, vi_it end, int rec_depth, int max_rec_depth) {
    int sz = end-begin;
    if(rec_depth > max_rec_depth) {
        if(sz > 16) base0_non_rec_heap_sort(begin, end);
        return;
    }

    if(sz <= 16) {
        insertion_sort(begin, end);
        return;
    }

    vi_it mid = partition(begin, end);
    intro_sort_rec(begin, mid, rec_depth+1, max_rec_depth);
    intro_sort_rec(mid+1, end, rec_depth+1, max_rec_depth);
}

void intro_sort(vi_it begin, vi_it end) {
    int sz = end-begin;
    if(sz <= 1) return;

    intro_sort_rec(begin, end, 1, ceil(log2(sz)));
    insertion_sort(begin, end);
}

int main() {
    Experiment experiment;
    Exp_Result res;
    vector<pair<sort_fcn, string>> nlogn_test_targets;
    vector<pair<sort_fcn, string>> n2_test_targets;

    nlogn_test_targets.push_back({merge_sort, "Merge Sort"});
    nlogn_test_targets.push_back({non_rec_merge_sort, "Non-recursive Merge Sort"});
    nlogn_test_targets.push_back({three_way_merge_sort, "3-way Merge Sort"});
    nlogn_test_targets.push_back({non_rec_three_way_merge_sort, "Non-recursive 3-way Merge Sort"});
    nlogn_test_targets.push_back({four_way_merge_sort, "4-way Merge Sort"});
    nlogn_test_targets.push_back({five_way_merge_sort, "5-way Merge Sort"});
    nlogn_test_targets.push_back({six_way_merge_sort, "6-way Merge Sort"});
    nlogn_test_targets.push_back({base0_heap_sort, "0-based Heap Sort"});
    nlogn_test_targets.push_back({base0_non_rec_heap_sort, "0-based Non-recursive Heap Sort"});
    nlogn_test_targets.push_back({base1_heap_sort, "1-based Heap Sort"});
    nlogn_test_targets.push_back({base1_non_rec_heap_sort, "1-based Non-recursive Heap Sort"});
    n2_test_targets.push_back({bubble_sort, "Bubble Sort"});
    n2_test_targets.push_back({insertion_sort, "Insertion Sort"});
    n2_test_targets.push_back({binary_insertion_sort, "Binary Insertion Sort"});
    n2_test_targets.push_back({selection_sort, "Selection Sort"});
    nlogn_test_targets.push_back({quick_sort, "Quick Sort"});
    n2_test_targets.push_back({library_sort, "Library Sort"});
    nlogn_test_targets.push_back({tim_sort, "Tim Sort"});
    nlogn_test_targets.push_back({no_galloping_tim_sort, "No-galloping Tim Sort"});
    n2_test_targets.push_back({cocktail_shaker_sort, "Cocktail Shaker Sort"});
    n2_test_targets.push_back({comb_sort, "Comb Sort"});
    nlogn_test_targets.push_back({tournament_sort, "Tournament Sort"});
    nlogn_test_targets.push_back({intro_sort, "Intro Sort"});

    int nlogn_sz = nlogn_test_targets.size();
    for(int i = 0; i < nlogn_sz; ++i) {
        res = experiment.nlogn_test(nlogn_test_targets[i].first);
        experiment.write_file(res, nlogn_test_targets[i].second);
        if(i < nlogn_sz-1) experiment.out_file << "\n\n" << endl;
    }

    int n2_sz = n2_test_targets.size();
    for(int i = 0; i < n2_sz; ++i) {
        res = experiment.n2_test(n2_test_targets[i].first);
        experiment.write_file(res, n2_test_targets[i].second);
        if(i < n2_sz-1) experiment.out_file << "\n\n" << endl;
    }

    experiment.out_file.close();
    return 0;
}