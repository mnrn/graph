/**
 * @brief Primのアルゴリズムのテストを行います
 * @note  関連URL: http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ALDS1_12_A
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <functional>
#include <utility>
#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>
#include <queue>
#include <iostream>



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



//****************************************
// 関数の定義
//****************************************

/**< @brief 優先度付きキューをmin優先度付きキューにするための < 演算子定義 */
static bool operator < (const edge& e, const edge& f) { return e.w > f.w; }



/**
 * @brief  Primのアルゴリズム
 *
 * @note   Primのアルゴリズムは、グラフの最短路を求めるDijkstraのアルゴリズムとほとんど同じように動作する
 *         Primのアルゴリズムは集合Aの辺が常に1つの木を形成するという性質を持つ. この木は任意の根rから開始し、Vの頂点全体を張るまで成長する
 *         各ステップでは、Aの頂点とある孤立点(Aの辺と接続していない頂点)を連結する軽い辺を木Aに加える
 *         Aに対して安全な辺だけがこの規則によってAに加えられるから、アルゴリズムが終了したとき、Aの辺は最小全域木を形成する
 *         各ステップでは木の重みの増加を限りなく小さく抑える辺を用いて木を成長させるので、これは貪欲戦略である
 *
 * @note   優先度付きキューの優先度更新を行わないため、優先度付きキューが空になるまでに行われる挿入の回数はΟ(E)であるが、
 *         EXTRACT-MIN呼び出し時に、黒頂点であれば無視をすることで、全体としての実行時間をΟ(ElgV)としている
 * 
 * @param  const graph_t& G グラフG
 * @param  index_t        r 最小全域木の根
 */
std::pair<edges_t, weight_t> prim(const graph_t& G, index_t r)
{
    std::int32_t n = G.size();
    std::size_t m = 0;
    for (auto&& es : G) { m += es.size(); }
    std::vector<std::int32_t> visited(n);
    edges_t A; weight_t w;

    
    for (index_t u = 0; u < n; u++) {
        visited[u] = false;  // 各頂点を白色に初期化
    }
    std::priority_queue<edge> Q;
    Q.push(edge(limits::nil, r, w = 0));    // 根rはキーを0に設定する
    while (!Q.empty()) {
        edge e = Q.top(); Q.pop();          // 軽い辺を取り出す
        index_t u = e.dst;
        if (visited[u]) { continue; }       // 取り出した辺が安全な辺ではない場合、再びループに戻り条件判定を行う
        
        for (auto&& f : G[u]) {             // uと隣接し、木に属さない各頂点vの更新を行う
            if (!visited[f.dst/*頂点v*/]) { Q.push(f); }
        }
        visited[u] = true;                  // 頂点uを黒色に彩色し、
        w += e.w;                           // 最小重みを更新する
        if (e.src != limits::nil) {  // アルゴリズムが終了したとき、min優先度付きキューは空であり、
            A.emplace_back(e);       // Gに対する最小全域木AはA = { (v, v.π) : v ∈ V - { r } }である  
        }
    }
    return std::make_pair(A, w);
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
    int n;   cin >> n;
    graph::graph_t G; G.resize(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k;
            cin >> k;
            if (k != -1) {
                G[i].emplace_back(i, j, k);
            }
        }
    }

    auto mst = prim(G, 0);
    cout << mst.second << endl;
    
    return 0;
}
