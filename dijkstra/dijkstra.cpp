/**
 * @brief  最短経路問題(shortest paths problem)における
 *         Dijkstraのアルゴリズム(Dijkstra's algorithm)を扱う
 * @date   2016/02/20 ~ 2016/03/13
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <queue>
#include <iostream>
#include "../graph/greedy.hpp"
#include "../graph/relax.hpp"
#include "dijkstra.hpp"



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 構造体の定義
//****************************************

struct state {
    index_t  u;  /**< G.Vに属する頂点u */
    weight_t d;  /**< 始点sからの距離d */

    /**< @brief <演算子オーバーロード */
    bool operator < (const state& s) const { return d > s.d; } // NOTE : min優先度付きキューのためにREVERSE

    state() = default;
    state(index_t u, weight_t d) : u(u), d(d) {}
};



//****************************************
// 関数の定義
//****************************************

/**
 * @brief  すべての辺重みが非負であるという仮定の下で、Dijkstra(ダイクストラ)のアルゴリズム(Dijkstra's algorithm)は
 *         重み付き有向グラフG = (V, E)上の単一始点最短路問題を解く. ここでは各辺(u, v) ∈ Eについてw(u, v) >= 0を仮定する
 *
 * @note   Dijkstraのアルゴリズムは、始点sからの最短路重みが最終的に決定された頂点の集合Sを管理する
 *         アルゴリズムは繰り返し、最小の最短路推定値を持つ頂点u ∈ V - Sを選択し、uをSに追加し、
 *         uから出るすべての辺を緩和する. ここではd値をキーとする頂点のmin優先度付きキューQを用いる
 *
 * @note   優先度付きキューの優先度更新を行わないため、優先度付きキューが空になるまでに行われる挿入の数はΟ(E)であるが、
 *         EXTRACT-MIN呼び出し時に、最短路の更新が行われないならば、無視をすることで、全体としての実行時間をΟ(ElgV)としている
 *
 * @param  const graph_t& G    非負の重み付き有向グラフG
 * @param  index_t        s    始点s
 * @return 始点sからの最短路重みが最終的に決定された頂点の集合S
 */
vertices_t dijkstra(const graph_t& G, index_t s)
{
    index_t n = G.size();
    vertices_t S(n);
    std::priority_queue<state> Q;

    
    initialize_single_source_with_color(S, s);     // すべての頂点のd値とπ値を初期化する
    Q.emplace(s, S[s].d);                          // このループの最初の実行ではu = sである
    while (!Q.empty()) {
        state p = Q.top(); Q.pop();
        index_t u = p.u; weight_t d = p.d;
        if (S[u].d < d) { continue; }
        for (auto&& e : G[u]) {        // 頂点uからでる辺(u, v)をそれぞれ緩和し、
            relax_with_heap(S, e, Q);  // uを経由することでvへの最短路が改善できる場合には、推定値v.dと先行点v.piを更新する
        }
        S[u].color = vcolor::black;    // 黒頂点は集合Sに属す
    }
    // 終了時点ではQ = φである. S = Vなので、すべての頂点u ∈ Vに対してu.d = δ(s, u)である
    // また、このとき、先行点部分グラフGπはsを根とする最短路木である
    return S;
}


/**
 * @brief  すべての辺重みが非負であるという仮定の下で、Dijkstra(ダイクストラ)のアルゴリズム(Dijkstra's algorithm)は
 *         重み付き有向グラフG = (V, E)上の単一始点最短路問題を解く. ここでは各辺(u, v) ∈ Eについてw(u, v) >= 0を仮定する
 *
 * @note   頂点v ∈ Vに対し、それぞれΟ(V)時間の操作を行うので、全体でΟ(V^2)時間を要す
 *
 * @param  const matrix_t& W    非負の重み付き有向グラフW
 * @param  index_t         s    始点s
 * @return 始点sからの最短路重みが最終的に決定された頂点の集合S
 */
vertices_t dijkstra(const matrix_t& W, index_t s)
{
    index_t n = W.size();
    vertices_t S(n);
    auto relax_pred = [](const vertices_t& V, index_t u) -> bool { return !V[u].visited; };    


    initialize_single_source_with_visitor(S, s);  // すべての頂点のd値とπ値を初期化する
    while(true) {     
        // 始点sからの最小の最短路推定値を持つ頂点u ∈ V - Sを選択する
        index_t u = extract_min(S, n);   
        if (u == limits::nil) { break; }    // 頂点uがNILを指すならば、探索は終了である
        for (index_t v = 0; v < n; v++) {
            relax(S, W, u, v, relax_pred);  // uを経由することでvへの最短路が改善できる場合には、推定値v.dと先行点v.piを更新する
        }
        S[u].visited = true;  // 黒頂点は集合Sに属す
    }
    // 終了時点では、S = Vなので、すべての頂点u ∈ Vに対してu.d = δ(s, u)である
    // また、このとき、先行点部分グラフGπはsを根とする最短路木である
    return S;
}



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END




//****************************************
// おまけ
//****************************************

/**
 * @brief  Dijkstraのアルゴリズムの詳細
 *
 * @note   Dijkstraのアルゴリズムは始点sからの最短路重みが最終的に決定された頂点の集合Sを管理する
 *         アルゴリズムは、繰り返し、最小の最短路推定値を持つ頂点u ∈ V - Sを選択し、uをSに追加し、uから出るすべての辺を緩和する
 *         以下に述べる実現はd値をキーとする頂点のmin優先度付きキューQを用いる
 *
 *         DIJKSTRA(G, w, s)
 *         1  INITIALIZE-SINGLE-SOURCE(G, s)
 *         2  S = φ
 *         3  Q = G.V
 *         4  while Q != φ
 *         5     u = EXTRACT-MIN(Q)
 *         6     S = S ∪ { u }
 *         7     for each vertex v ∈ G.Adj[u]
 *         8         RELAX(u, v, w)
 *
 *         第1行でdとπ値をいつもどおりに初期化し、第2行で集合Sを空集合に初期化する
 *         第4~8行のwhile文の各繰り返しの開始直前では、アルゴリズムは不変式Q = V - Sを維持する
 *         第3行ではmin優先度付きキューQをVのすべての頂点を含むように初期化する. この時点ではS = φだから、第3行を実行した後では不変式Q = V - Sから抽出し、
 *         集合Sに挿入するから、不変式を維持できる.(このループの最初の実行ではu = sである.) 頂点uはV - Sに属する頂点の中で最小の最短路推定値を持つ. 
 *         つぎに、第7~8行ではuから出る辺(u, v)のそれぞれを緩和し、uを経由することでvへの最短路が改善できる場合には、推定値v.dと先行点v.πを更新する
 *         どの頂点も第3行以降でQに挿入されることはなく、各頂点はちょうど1回だけQから抽出されてSに挿入されるので、第4~8行のwhileはちょうど|V|回だけ
 *         繰り返されることを確認せよ
 *         
 * @note   DijkstraのアルゴリズムはV - Sの中で常に"最も軽い"あるいは"最も近い"頂点を集合Sに挿入するから、貪欲戦略に基づいている
 *         貪欲戦略は一般には最適解を保証しないが、Dijkstraのアルゴリズムは実際に最短路を計算する
 *         ある頂点uを集合Sに挿入するときには常にu.d = δ(s, u)であることを示すことがキーになる
 *
 * @note   Dijktraのアルゴリズムは、幅優先探索と最小全域木を求めるPrimのアルゴリズムの両方と類似点を持っている
 *         集合Sと幅優先探索の黒頂点集合との対応という点でDijkstraのアルゴリズムは幅優先探索と似ている
 *         Sの頂点が最終的な最短路重みを持つように、幅優先探索の黒頂点も正しい幅優先距離をもつ
 *         DijkstraのアルゴリズムとPrimのアルゴリズムの類似点は、ともに、min優先度付きキューを用いて与えられた集合(Dijkstraのアルゴリズムでは集合S,
 *         Primのアルゴリズムでは成長中の木)に属さない"最も軽い頂点"を求め、この頂点を集合に加え、この集合に属さない頂点の重みを適切に調節するところにある
 */




/**
 * @brief    おまけ...boostのサイトに存在するDijkstraのアルゴリズムに対する疑似コード
 *
 * @details  DIJKSTRA(G, s, w)
 *           1  for each vertex u in V
 *           2      d[u] := infinity
 *           3      p[u] := u
 *           4      color[u] := WHITE
 *           5  end for
 *           6  color[s] := GRAY
 *           7  d[s] := 0
 *           8  INSERT(Q, s)
 *           9  while (Q != φ)
 *           10     u := EXTRACT-MIN(Q)
 *           11     S := S ∪ { u }
 *           12     for each vertex v in Adj[u]
 *           13         if (w(u, v) + d[u] < d[v])
 *           14             d[v] := w(u, v) + d[u]
 *           15             p[v] := u
 *           16             if (color[v] = WHITE)
 *           17                 color[v] := GRAY
 *           18                 INSERT(Q, v)
 *           19             else if (color[v] = GRAY)
 *           20                 DECREASE-KEY(Q, v)
 *           21             else
 *           22                 ...
 *           23     end for
 *           24     color[u] := BLACK
 *           25 end while
 *           26 return (d, p)
 *
 * @note     同サイトのPrimのアルゴリズムの疑似コードと同様、重要な部分を含んでいる
 *           やはり問題は、20行目のDECREASE-KEY呼び出しであるが、これもPrimのアルゴリズムの場合と同様の理由で
 *           16~20行目の操作を
 *           ex1  if (color[v] != BLACK)
 *           ex2    INSERT(Q, v)
 *           として構わない. この場合、優先度付きキューが空になるまでに挿入される頂点数はΟ(E)であると考えられるため、
 *           第9行のループ回数もΟ(E)であり、DIJKSTRAのアルゴリズムの総実行時間はΟ(ElgE)となることがわかる
 *           |E| < |V|^2に注意するとDIJKSTRAの総実行時間をΟ(ElgV)と書き直すことができる
 */






















