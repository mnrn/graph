/**
 * @brief  最大フローを求めるFord-Fulkersonアルゴリズムを扱う
 *
 * @note   フローネットワーク(flow network)G = (V, E)は、各辺(u, v) ∈ Eが非負の容量(capacity)c(u, v) >= 0を持つ有向グラフである
 *         さらに、Eが辺(u, v)を含むならば、逆向きの辺(v, u)を含まない.
 *         (u, v) !∈ Eならば、便宜上、c(u, v) = 0と定義し、自己ループは許さない
 *         フローネットワークには特別な2つの頂点、入口(source)sと出口(sink)tが指定されている. 便宜上、各頂点は入口から出口へのある道上にあると仮定する
 *         すなわち、各頂点v ∈ Vに対して、フローネットワークはある道s~>v~>tを含む. したがって、グラフは連結であり、s以外のすべての頂点には少なくとも
 *         1本の辺が入ってくるので、|E| >= |V| - 1である
 *
 *         G = (V, E)を容量関数cを持つフローネットワークとする. sをこのネットワークの入口, tを出口とする
 *         Gにおけるフロー(flow)はつぎの2条件を満たす実数値関数f: V x V -> Rである
 *           容量制限(capacity constraint)  : すべてのu,v ∈ Vに対して、0 <= f(u, v) <= c(u, v)でなければならない
 *           フロー保存則(flow conservation) : すべてのu ∈ V - { s, t }に対して、Σ[v∈V]f(v, u) = Σ[v∈V]f(u, v)でなければならない
 *                                           (u, v) !∈ Eならば、uからvへのフローは存在せず、f(u, v) = 0である
 *
 *         非負の値を取る量f(u, v)を頂点uから頂点vへのフロー(flow)と呼ぶ. フローfの値(value)|f|を
 *           |f| = Σ[v∈V]f(s, v) - Σ[v∈V]f(v, s)、
 *         すなわち、入口から流れ出るフローの合計と入口に流れ込むフローの合計の差として定義する(ここで、記号|・|でフローの値を表すが、これは絶対値でも要素数でもない)
 *         フローネットワークは一般に入口に入る辺を持たず、和Σ[v∈V]f(v, s)が与える入口へ流入するフローは0である
 *         しかし、残余ネットワークでは入口へ流入するフローが重要になるので、これを含めて定義する
 *         最大フロー問題(maximum-flow problem)は、入口sと出口tを持つフローネットワークGが与えられたとき、sからtへの最大の値をもつフローを求める問題である
 *
 *         フローの定義の2つの条件を説明しておく. 容量制限は、ある頂点から別の頂点へのフローは非負でなければならず、また与えられた容量を超えてはいけないという制約である
 *         フロー保存則は、入口と出口以外の頂点では、流れ込むフローの合計が流れ出るフローの合計に一致すること、すなわち「流入量 = 流出量」が成立することを要請する
 *
 * @date   2016/02/22 ~ 2017/01/03
 */



//****************************************
// インクルードガード
//****************************************

#ifndef __FORD_FULKERSON_HPP__
#define __FORD_FULKERSON_HPP__



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include "../graph/graph.hpp"



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN


//****************************************
// 構造体の定義
//****************************************

/**
 * @brief   基本Ford-Fullkersonアルゴリズムを扱う
 *
 * @details FORD-FULKERSON-METHOD(G, s, t)
 *          1  フローfを0に初期化する
 *          2  while 増加可能経路pが残余ネットワークGfに存在する
 *          3    フローfをpに沿って増やす
 *          4  return f
 *
 *          上で与えたFORD-FULKERSON-METHODの疑似コードを展開したのが、以下のFORD-FULKERSONアルゴリズムである
 *
 *          FORD-FULKERSON(G, s, t)
 *          1  for 各辺(u, v) ∈ G.E
 *          2      (u, v).f = 0
 *          3      (v, u).f = 0
 *          4  while 残余ネットワークGfにsからtへの道pが存在する
 *          5      cf(p) = min{cf(u, v) : (u, v)はpに属する }
 *          6      for 各辺(u, v) in p
 *          7          (u, v).f = (u, v).f + cf(p)
 *          8          (v, u).f = (v, u).f - cf(p)
 *
 *          第1~3行目でフローfを0に初期化する. 第4~8行のwhile文ではGf上の増加可能経路pを見つけ、pに沿ってフローfを残余容量cf(p)だけ増やす操作を繰り返す
 *          道p上の各残余辺は元のネットワークの辺か、その逆向き辺である. 第6~8行では適切にフローを更新する. 増加可能経路がなければ、フローfは最大フローである
 */
struct ford_fulkerson {
    stamps_t visited;             /**< すでに訪問済みか？ */
    matrix_t c, f;                /**< 辺(u, v) ∈ Eの容量属性(u, v).cとフロー属性(u, v).f */
    std::vector<indices_t> Gf;    /**< 残余ネットワークGf */
    index_t n;                    /**< 頂点v ∈ Vの数 */
    capacity_t augment;           /**< フローの増加数 */

    explicit ford_fulkerson(std::size_t size) : c(size, array_t(size)), f(size, array_t(size)),
                                   Gf(size, indices_t()), n(size), augment(0) {}
    explicit ford_fulkerson(const graph_t& G)
    {
        *this = ford_fulkerson(G.size());
        for (index_t i = 0; i < n; ++i) {
            for (auto&& e : G[i]) { add_edge(e.src, e.dst, e.c); }
        }
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
     * @brief  Ford-Fulkersonのアルゴリズムを実行する
     * @note   アルゴリズムの実行時間はΟ(E|f*|)であるが、コンストラクタで隣接行列を生成しているので全体の実行時間はΟ(E|*f| + V^2)
     * @param  フローネットワークの入口(source) s
     * @param  フローネットワークの出口(sink)   t
     * @return フローネットワークの最大フロー
     */
    capacity_t compute(index_t s, index_t t)
    {
        capacity_t flow = 0;    // execute終了時にsからtへの最大フローとなる値

        // NOTE : コンストラクタ呼び出し時にFORD-FULKERSONの第1~3行の初期化と同様の操作は終了している
        
        // NOTE : 増加可能経路pをDFSで辿りつつ、FORD-FULKERSONの第5~8行を実行する
        while (dfs(s, t)) {
            flow += augment;
        }
        return flow;
    }

    /**
     * @brief  深さ優先探索を用いて残余ネットワークGfにsからtへの道(増加可能経路(augment path))pを探索し、
     *         増加可能経路pに沿ってフローfを残余容量cf(p)だけ増やす
     * @param  index_t s 残余ネットワークGfの頂点u
     * @param  index_t s フローネットワークの出口(sink) t
     * @param  capacity_t flow 入口sから現在探索している頂点uまで道qの残余容量cf(q)
     * @return capacity_t cf_p
     */
    capacity_t dfs_visit(index_t u, index_t t, capacity_t flow)
    {
        visited[u] = true;           // 訪問印を刻む
        if(u == t) { return flow; }  // 出口(sink)tに達した場合、再帰は底をつく
        
        for(auto&& v : Gf[u]) {  // 各頂点v ∈ Adj[u]を吟味するので、深さ優先探索は辺(u, v)を探索する(explore)という
            // vが白ではない、または残余容量がゼロならば辺(u, v)を調べる必要はない
            if(visited[v] || cf(u, v) == 0) { continue; }
            
            // 再帰的にDFS-VISITを呼び出し、残余容量cf(p)を得る
            capacity_t cf_p = dfs_visit(v, t, std::min(flow, cf(u, v)));
            
            if (cf_p > 0) {  // cf(p)がゼロでないならば、
                f[u][v] += cf_p;  // 元のネットワークの辺のフローを加え、
                f[v][u] -= cf_p;  // 逆向き辺のフローを引く
                return cf_p;
            }
        }
        return 0;
    }

    /**
     * @brief  深さ優先探索を用いて残余ネットワークGfにsからtへの道(増加可能経路(augment path))pを探索し、
     *         増加可能経路pに沿ってフローfを残余容量cf(p)だけ増やす
     * @param  index_t s 残余ネットワークGfの頂点u
     * @param  index_t s フローネットワークの出口(sink) t
     * @return bool_t  残余ネットワークGfにsからtへの道(増加可能経路(augment path))pが存在するか？
     */
    bool_t dfs(index_t u, index_t t)
    {
        visited.assign(n, false);  // DFS-VISITでpを探す前に、訪問フラグをすべてfalseにする
        augment = dfs_visit(u, t, limits::inf);
        return augment > 0;        // pが存在するか？
    }

    /**< @brief 頂点対u, v ∈ Vにおける残余容量cf(u, v)を返す */
    capacity_t cf(index_t u, index_t v) const
    {
        return c[u][v] - f[u][v];
    }
};



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



#endif  // end of __FORD_FULKERSON_HPP__

