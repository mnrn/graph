/**
 * @brief 幅優先探索のテスト
 * @note  関連URL: http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ALDS1_11_C
 * @date  2016/02/13
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <iostream>
#include <iterator>
#include <functional>
#include <utility>
#include <algorithm>
#include <vector>
#include <cassert>
#include <queue>



//****************************************
// オブジェクト形式マクロの定義
//****************************************

#define GRAPH_BEGIN namespace graph {
#define GRAPH_END   }



//****************************************
// 構造体の定義
//****************************************

GRAPH_BEGIN

using weight_t   = std::int32_t;  /**< 辺(u, v)への重みwを表す型 */
using index_t    = std::int32_t;  /**< 頂点vの添字を表す型       */
using capacity_t = weight_t;      /**< 辺(u, v)の容量を表す型    */



namespace limits {
    enum {  // namespace graphによってscopedではあるが強く型付けされたenum(strongly-typed enum)ではない
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

    bool operator < (const edge& e) const { return w < e.w; }
    bool operator > (const edge& e) const { return w > e.w; }
};



using edges_t    = std::vector<edge>;      /**< グラフG=(V, E)の辺集合E   */
using vertices_t = std::vector<vertex>;    /**< グラフG=(V, E)の頂点集合V */
using array_t    = std::vector<weight_t>;  /**< 重みwの配列  */
using indices_t  = std::vector<index_t>;   /**< 頂点の添字配列 */
using matrix_t   = std::vector<array_t>;   /**< グラフGの隣接行列表現(および表行列表現) */
using graph_t    = std::vector<edges_t>;   /**< グラフGの隣接リスト表現(こちらを主に使用する) */



//****************************************
// 関数の定義
//****************************************

/**
 * @brief  幅優先探索を行います
 *
 * @note   幅優先探索(breadth-first-search)は最も単純なグラフ探索アルゴリズムの1つである
 *         グラフG=(V, E)を始点(source vertex)sが与えられたとき、幅優先探索はGの辺を組織的に探索して、
 *         sから到達可能なすべての頂点を"発見し"、すべての到達可能な頂点についてsからの距離(変数の最小値)を計算する
 *         さらに、sを根とし、到達可能な頂点をすべて含む「幅優先木」を構成する. sから到達可能な任意の頂点vについて、
 *         幅優先木におけるsからvへの単純道はGにおけるsからvへの「最短路」、すなわち最小数の辺を含む道に対応する
 *         有向グラフと無向グラフのどちらに対してもこのアルゴリズムは正しく動く
 *
 *         探索頂点と未探索頂点の境界を、境界の幅一杯にわたって一様に拡張することが幅優先探索という名前の由来である
 *         すなわち、このアルゴリズムはsから距離k+1にある頂点を発見する前に距離kにあるすべての頂点を発見する
 *
 * @note   手続きBFSはグラフを探索しながら幅優先木を構築する.この木はπ属性に対応する.形式的に言うと、
 *         sを始点とするグラフG=(V, E)に対して、Gの先行点部分グラフ(predecessor subgraph)をGπ=(Vπ, Eπ)として定義する.ただし、
 *           Vπ = { v ∈ V : v.π != NIL } ∪ {s}
 *         かつ
 *           Eπ = { (v.π, v) : v ∈ Vπ - {s} }
 *         である.Vπがsから到達可能な全頂点から構成され、すべてのv∈Vπに対して、sからvに至る唯一の単純道がGπに存在し、
 *         しかもこれがGにおけるsからvに至る最短路になっているとき、先行点部分グラフGπを幅優先木(breadth-first tree)と呼ぶ
 *         幅優先木は連結で|Eπ| = |Vπ| - 1を満たすから、実際に木である. Eπの辺を木辺(tree edge)と呼ぶ
 *
 * @note   BFSの総実行時間はΟ(V+E)である.したがって、幅優先探索はGの隣接リスト表現のサイズの線形時間で走る
 *
 * @param  const graph& G  グラフG
 * @param  index_t s  始点s
 * @return 幅優先木
 */
vertices_t bfs(const graph_t& G, index_t s)
{
    index_t n = G.size();
    vertices_t V(n);

    for (auto&& u : V) {            // すべての頂点uについて、
        u.color = vcolor::white;    // uを白に彩色し、
        u.d     = limits::inf;      // u.dを無限大に設定し、
        u.pi    = limits::nil;      // uの親をNILに設定する
    }
    // 手続き開始と同時に始点sを発見すると考え、
    V[s].color = vcolor::gray;      // 始点sを灰色に彩色する 
    V[s].d     = 0;                 // s.dを0に初期化し、
    V[s].pi    = limits::nil;       // 始点の先行点をNILに設定する

    std::queue<index_t> Q;
    Q.push(s);                      // sだけを含むようにQを初期化する

    // 以下のwhile文に対して、つぎのループ不変式が成立する
    // while文の条件判定を行う時点ではキューQはすべての灰頂点を含む
    while (!Q.empty()) {
        index_t u = Q.front(); Q.pop();
        for (auto&& e : G[u]) {                  // uの隣接リストに
            index_t v = e.dst;                   // 属する各頂点vを考える
            if (V[v].color == vcolor::white) {   // vが白ならvは未発見である
                V[v].color = vcolor::gray;       // vを灰色に彩色し、
                V[v].d     = V[u].d + 1;         // 距離v.dをu.d+1に設定し、
                V[v].pi    = u;                  // uをその親v.piとして記録し、
                Q.push(v);                       // vをキューQの末尾に置く
            }
        }
        V[u].color = vcolor::black;  // uの隣接リストに属するすべての頂点の探索が完了すると、この頂点を黒に彩色する
    }
    // ある頂点を灰に彩色したときには、この頂点をQへ挿入し、ある頂点をQから削除したときには、この頂点を黒に彩色するので、
    // ループ不変式が保存される

    return V;
}


GRAPH_END



int main()
{
    using namespace graph;
    graph_t G;
    int u, k, v, n;
    std::cin >> n;
    G.resize(n);
    for (int i = 0; i < n; i++) {
        std::cin >> u >> k;
        u--;
        G[i].resize(k);
        for (int j = 0; j < k; j++) {
            std::cin >> v;
            v--;
            G[i].emplace_back(i, v);
        }
    }

    auto bfstree = bfs(G, 0);
    for (int i = 0; i < n; i++) {
        std::cout << i + 1 << " "
                  << ((bfstree[i].d == limits::inf)
                      ? -1 : bfstree[i].d) << std::endl;
    }
    return 0;
}

