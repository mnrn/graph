/**
 * @brief  最大フローを求めるEdmondsKarpアルゴリズムを扱う
 * @date   2016/03/03 ~ 2017/01/03
 */



//****************************************
// インクルードガード
//****************************************

#ifndef __EDMONDS_KARP_HPP__
#define __EDMONDS_KARP_HPP__



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include "../graph/graph.hpp"
#include <queue>



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 構造体の定義
//****************************************

/**
 * @brief Edmonds-Karpのアルゴリズム
 * @note  Ford-Fulkerson法において、増加可能経路を幅優先探索を用いて探索することでFORD-FULKERSONの計算時間の上界を改善できる
 *        すなわち、残余ネットワークの中でsとtを結ぶ最短路を増加可能経路として選択するのである. ただし、残余ネットワークの各辺(u, v)の距離(重み)は1である
 *        Ford-Fulkerson法をこのように実現したものをEdmonds-Karpアルゴリズム(Edmonds-Karp algorithm)と呼ぶ
 */
struct edmonds_karp {
    indices_t pi;                 /**< 頂点vの先行点属性 */
    stamps_t visited;             /**< すでに訪問済みか？ */
    matrix_t c, f;                /**< 辺(u, v) ∈ Eの容量属性(u, v).cとフロー属性(u, v).f */
    std::vector<indices_t> Gf;    /**< 残余ネットワークGf */
    index_t n;                    /**< 頂点v ∈ Vの数 */

    explicit edmonds_karp(std::size_t size) :
        pi(size), visited(size), c(size, array_t(size)), f(size, array_t(size)),
        Gf(size, indices_t()), n(size) { }
    explicit edmonds_karp(const graph_t& G)
    {
        *this = edmonds_karp(G.size());
        for (index_t i = 0; i < n; ++i) { for (auto&& e : G[i]) { add_edge(e.src, e.dst, e.c); } }
    }

    /**
     * @brief  容量cおよびフローfの初期化と残余ネットワークGfの生成部分
     * @note   入口sと出口tを持つフローネットワークをG = (V, E)とする. fをGのフローとし、頂点対u, v ∈ Vを考える
     *         (u, v)の残余容量(residul capacity)cf(u, v)を
     *           cf(u, v) = { c(u, v) - f(u, v)  (u, v) ∈ Eのとき
     *                        f(v, u)            (v, u) ∈ Eのとき
     *                        0                  それ以外         }
     *         と定義する. 定義から(u, v) ∈ Eならば、(v, u) !∈ Eだから、任意の頂点対に対して上記の式のちょうど1つの場合が対応する
     *
     *         フローネットワークG = (V, E)とフローfが与えられたとき、fによって誘導される残余ネットワーク(residual network)は
     *           Ef = { (u, v) ∈ V x V : cf(u, v) > 0 }
     *         によって定義されるGf = (V, Ef)である
     *         残余ネットワークの各辺、すなわち残余辺(residual edge)には正のフローを流すことができる
     *         Efの辺はEの辺かその逆向きの辺であり、したがって |Ef| <= 2|E|である
     */
    void add_edge(index_t u, index_t v, capacity_t cap)
    {
        c[u][v] = cap; f[u][v] = 0;              // cおよびfの初期化
        c[v][u] = 0;   f[v][u] = 0;
        Gf[u].push_back(v); Gf[v].push_back(u);  // Gfの生成
    }


    /**
     * @brief  Edmonds-Karpのアルゴリズムを実行する
     * @note   アルゴリズムの実行時間はΟ(VE^2)であるが、コンストラクタで隣接行列を生成しているので全体の実行時間はΟ(VE^2 + V^2)
     * @param  index_t s フローネットワークの入口(source) s
     * @param  index_t t フローネットワークの出口(sink) t
     * @return フローネットワークの最大フロー
     */
    capacity_t compute(index_t s, index_t t)
    {
        capacity_t flow = 0;  // execute終了時にsからtへの最大フローとなる値

        // NOTE: コンストラクタ呼び出し時にFORD-FULKERSONの第1~3行の初期化と同様の操作は終了している
        
        while (bfs(s, t)) { flow += proc(s, t); }  // BFSでpを探し、pが存在したならば、フローを更新する
        return flow;
    }


    /**
     * @brief  幅優先探索を用いて残余ネットワークGfにsからtへの道(増加可能経路(augment path))pを探索する
     * @param  index_t s フローネットワークの入口(source) s
     * @param  index_t s フローネットワークの出口(sink) t
     * @return 増加可能経路pが存在するか否か
     */
    bool_t bfs(index_t s, index_t t)
    {
        pi.assign(n, limits::nil);
        visited.assign(n, false);
        std::queue<index_t> Q;

        // 手続き開始と同時に始点sを発見したと考え、
        pi[s] = limits::nil;  // 始点の先行点をNILで初期化する
        visited[s] = true;    // 訪問印を刻む
        
        Q.push(s);  // sだけを含むようにキューを初期化する
        while (!Q.empty()) {
            index_t u = Q.front(); Q.pop();
            for (auto&& v : Gf[u]) {
                // vが白でない、または残余容量がゼロならば、辺(u, v)を調べる必要がない
                if (visited[v] || cf(u, v) == 0) { continue; }
                
                // 上記の条件にいずれも当てはまらない場合、
                visited[v] = true;  // 訪問印を刻み、
                pi[v] = u;          // uをその親v.piとして記録し、
                Q.push(v);          // vをキューQの末尾に置く
            }
        }
        return visited[t];  // 最終的な結果は出口節点tを訪問するか、pが空のどちらかである
    }


    /**
     * @brief  増加可能経路pに沿ってフローfを残余容量cf(p)だけ増やす
     * @param  index_t s フローネットワークの入口(source) s
     * @param  index_t t フローネットワークの出口(sink) t
     * @return capacity_t cf_p
     */
    capacity_t proc(index_t s, index_t t)
    {
        capacity_t cf_p = limits::inf;
        index_t v = t;

        // フローを増やす量を決定する. sinkからsourceへと計算した経路で最小のものと等しい
        while (v != s) {
            cf_p = std::min(cf_p, cf(pi[v], v));
            v = pi[v];
        }

        // 経路上で得られた最小値分だけ増加させる
        v = t;
        while (v != s) {
            f[pi[v]][v] += cf_p;  // 元のネットワークの辺(前方辺)のフローを加え、
            f[v][pi[v]] -= cf_p;  // 逆向き辺(後方辺)のフローを引く
            v = pi[v];
        }
        
        return cf_p;
    }

    /**< @brief 頂点対(u, v) ∈ Vにおける残余容量cf(u, v)を返す */
    capacity_t cf(index_t u, index_t v) const
    {
        return c[u][v] - f[u][v];
    }
};



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



#endif  // end of __EDMONDS_KARP_HPP__

