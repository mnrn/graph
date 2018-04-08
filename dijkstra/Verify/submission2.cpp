/**
 * @brief  最短経路問題(shortest paths problem)におけるDijkstraのアルゴリズム(Dijkstra's algorithm)を扱う
 *
 * @note   関連URLその1: http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ALDS1_12_B
 *               その2: http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ALDS1_12_C
 *
 * @date   2016/02/20
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>
#include <cassert>
#include <functional>
#include <queue>



//****************************************
// オブジェクト形式マクロの定義
//****************************************

#define GRAPH_BEGIN namespace graph {
#define GRAPH_END   }



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 構造体の定義
//****************************************

using weight_t   = std::int32_t;  /**< 辺(u, v)への重みwを表す型 */
using index_t    = std::int32_t;  /**< 頂点vの添字を表す型       */
using capacity_t = weight_t;      /**< 辺(u, v)の容量を表す型    */



namespace limits {
    enum {  // scopedではあるが強く型付けされたenum(strongly-typed enum)ではない
        inf = std::numeric_limits<weight_t>::max() / 3,  /**< @brief 辺が存在しない場合に使用される値     */
        nil = std::numeric_limits<index_t>::min() / 3,   /**< @brief 先行点が存在しない場合に使用される値 */
    };
}



/**
 * @brief  頂点色列挙構造体(scoped enum)
 * @detail (u, v) ∈ Eで頂点uが黒ならば頂点vは灰か黒である
 *         すなわち、黒頂点に隣接する全ての頂点は発見済みである
 *         灰頂点は白頂点に隣接することがあり、これらの頂点が既発見頂点と未発見頂点の境界をなす
 */
enum struct vcolor : std::int32_t {
    white,  /**< 未発見頂点 */
    black,  /**< 既発見頂点 */
    gray,   /**< 既見済頂点 */
};

/**
 * @brief グラフ用ノード(頂点) 
 */
struct vertex {
    union {
        weight_t d;       /**< 始点sからの距離  */
        weight_t key;     /**< Primのアルゴリズムにおいて木に属するある頂点とを結ぶ重み */
    };
    index_t pi;           /**< 先行頂点(の添字) */
    union {
        vcolor color;     /**< 頂点の色        */
        bool  visited;    /**< 発見済みか?     */
    };
    // weight_t f;           /**< 終了時刻印(DFSにおいて、黒色に彩色されたとき、刻まれる)     */
    vertex() : d(0), pi(0), color(vcolor::white)/*, f(0)*/ {}
};

/**
 * @brief グラフ用エッジ(辺)
 * @note  G = (V, E)を重み関数wを持つ重み付きグラフとすると、
 *        辺(u, v) ∈ Eの重みはw(u, v)と表される
 */
struct edge {
    index_t  src;   /**< 辺の始点u */
    index_t  dst;   /**< 辺の終点v */
    union {
        weight_t w;     /**< 辺(u, v)への重み(コスト) */
        capacity_t c;   /**< 辺(u, v)の容量 */
    };
    edge() = default;
    edge(index_t src, index_t dst) : src(src), dst(dst), w(1) {}
    edge(index_t src, index_t dst, weight_t w) : src(src), dst(dst), w(w) {}
};



using edges_t    = std::vector<edge>;      /**< グラフG=(V, E)の辺集合E   */
using vertices_t = std::vector<vertex>;    /**< グラフG=(V, E)の頂点集合V */
using array_t    = std::vector<weight_t>;  /**< 重みwの配列  */
using indices_t  = std::vector<index_t>;   /**< 頂点の添字配列 */
using matrix_t   = std::vector<array_t>;   /**< グラフGの隣接行列表現(および表行列表現) */
using graph_t    = std::vector<edges_t>;   /**< グラフGの隣接リスト表現(こちらを主に使用する) */



struct state {
    index_t  u;  /**< G.Vに属する頂点u */
    weight_t d;  /**< 始点sからの距離d */

    /**< @brief <演算子オーバーロード */
    bool operator < (const state& s) const { return d > s.d; } // NOTE : min優先度付きキューのためREVERSE

    state() = default;
    state(index_t u, weight_t d) : u(u), d(d) {}
};




//****************************************
// 関数の定義
//****************************************

/**
 * @brief  Θ(V)の手続きによって最短路推定値と先行点および頂点色を初期化する
 * @note   初期化の後、すべてのv ∈ Vについてv.π = NIL、
 *         すべてのv ∈ V - {s}についてv.d = ∞、 v.color = WHITEである
 */
static inline void initialize_single_source_with_color(vertices_t& V, index_t s)
{
    for (auto&& v : V) {
        v.d     = limits::inf;
        v.pi    = limits::nil;
        v.color = vcolor::white;
    }
    V[s].d = 0;
    V[s].color = vcolor::gray;
}


/**
 * @brief  辺(u, v)を緩和すると同時に、頂点vおよび道s~>vの重みをmin優先度付きキューQに挿入する
 * 
 * @tparam PriorityQueue min優先度付きキューの型
 * @param vertices_t&    V 頂点集合V
 * @param const edge&    e 辺(u, v)
 * @param PriorityQueue& Q min優先度付きキュー
 */
template<class PriorityQueue>
void relax_with_heap(vertices_t& V, const edge& e, PriorityQueue& Q)
{
    index_t u = e.src, v = e.dst;
    if (V[v].color != vcolor::black && V[v].d > V[u].d + e.w) {
        V[v].d     = V[u].d + e.w;
        V[v].pi    = u;
        V[v].color = vcolor::gray;
        Q.emplace(v, V[v].d);
    }
}


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



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



//****************************************
// エントリポイント
//****************************************

int main(void)
{
    using namespace std;
    using namespace graph;
    graph_t G;
    int n, k, u, v, c;
    
    cin >> n; G.resize(n);
    for (int i = 0; i < n; i++) {
        cin >> u >> k;
        G[u].resize(k);
        for (int j = 0; j < k; j++) {
            cin >> v >> c;
            G[u][j].src = u;
            G[u][j].dst = v;
            G[u][j].w   = c;
        }
    }

    vertices_t S = dijkstra(G, 0);
    for (int i = 0; i < n; i++) {
        cout << i << " " << S[i].d << endl; 
    }
    return 0;
}
