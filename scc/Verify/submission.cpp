/**
 * @brief 強連結(分解)アルゴリズムのテスト
 * @date  2016/03/12~2017/01/03
 * @note  関連URL: http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ALDS1_11_D
 */


//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <iterator>
#include <functional>
#include <utility>
#include <algorithm>
#include <vector>
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
 * @brief  トポロジカルソートを行います
 *
 * @note   有向非巡回グラフG = (V, E)のトポロジカルソート(topological sort)は頂点集合上の線形順序で、
 *         Gが辺(u, v)を含むならば、この線形順序でuがvより先に現れるものである. (グラフに巡回路があればこのような線形順序は存在しない)
 *         グラフのトポロジカルソートは、すべての有向辺が左から右へ向かう、水平線上での頂点の並べ方である
 *
 * @note   つぎの簡単なアルゴリズムは有向非巡回グラフをトポロジカルソートする
 *
 *         TOPOLOGICAL-SORT(G)
 *         1 各頂点の終了時刻v.fを計算するためにDFS(G)を呼び出す
 *         2 ある頂点の探索が終了するたびに、この頂点を連結リストの先頭に挿入する
 *         3 return 頂点の連結リスト
 *
 * @note   深さ優先探索にΘ(V + E)時間かかり、|V|個の頂点のそれぞれを連結リストの先頭に挿入するのにΟ(1)時間かかるので、
 *         トポロジカルソートはΘ(V + E)時間で実行できる
 *
 * @param  const graph_t& G 有向非巡回グラフ
 * @return 既ソートリスト
 */
array_t tsort(const graph_t& G)
{
    index_t n = G.size();
    std::vector<vcolor> color(n, vcolor::white);
    array_t lst(n);

    // 白節点を訪れる
    std::function<bool(index_t)> dfs_visit = [&](index_t u) {
        color[u] = vcolor::gray;    // uを灰に彩色する
        for (auto&& e : G[u]) {     // vと隣接する各頂点wを調べ、
            index_t w = e.dst;
            if (color[w] == vcolor::white
             && !dfs_visit(w)) { return false; }  // wが白なら再帰的にwを調べる
        }
        color[u] = vcolor::black;   // uを黒に彩色する
        lst.push_back(u);           // リストの末尾に挿入する
        return true;
    };


    // 各頂点vの終了時刻v.fを計算するためにDFS(G)を呼び出す
    for (index_t v = 0; v < n; ++v) {
        if (color[v] == vcolor::white && !dfs_visit(v)) { return {}; };
    }
    reverse(lst.begin(), lst.end());  // リストが逆順にソートされているのでreverseを行う
    return lst;     // 頂点のリストを返す
}




/**
 * @brief  強連結成分(分解)アルゴリズム
 *
 * @note   グラフを強連結成分に分解した後、個々の強連結成分上でアルゴリズムを実行し、
 *         得られた解を成分間の連結構造にしたがって組み合わせて最終的に解を得る
 *
 * @note   グラフG = (V, E)の強連結成分を求めるアルゴリズムはGの転置を用いる
 *         Gの転置はグラフG^T = (V, E^T), E^T = { (u, v) : (v, u) ∈ E }である
 *         すなわち、E^TはGの辺の方向を逆にしたものである. Gの隣接リストが与えられると、
 *         G^TをΟ(V + E)時間で生成できる. 興味深いことに、GとG^Tはまったく同じ強連結成分を持つ
 *         G上でuとvが互いに到達可能であることと、G^T上でこれらが互いに到達可能であることとは等価である
 *
 * @note   以下の線形時間(すなわち、Θ(V + E)時間)アルゴリズムは、深さ優先探索を2回、最初にG上で、
 *         つぎにG^T上で実行することによって有向グラフG = (V, E)の強連結成分を求める
 *
 *         STRONGLY-CONNECTED-COMPONENTS(G)
 *         1 DFS(G)を呼び出し、各頂点uに対して終了時刻u.fを計算する
 *         2 G^Tを計算する
 *         3 DFS(G^T)を呼び出すが、DFSの主ループでは(第1行で計算した)u.fの降順で頂点を探索する
 *         4 第3行で生成した深さ優先森の各木の頂点を、それぞれ分離された強連結成分として出力する
 *
 * @param  const graph_t& G グラフG
 * @return components[v] 頂点vが含まれる連結成分の番号となるような集合
 */
indices_t scc(const graph_t& G)
{
    index_t n = G.size(); vertices_t vs(n);
    indices_t components(n, -1);
    std::vector<vcolor> color(n, vcolor::white);
    graph_t GT(n);

    // 再帰的に白頂点を訪れる
    std::function<void(index_t, index_t)> dfs_visit = [&](index_t u, index_t k) {
        color[u] = vcolor::gray;
        components[u] = k;
        for (auto&& e : GT[u]) {
            index_t w = e.dst;
            if (color[w] == vcolor::white) {
                dfs_visit(w, k);
            }
        }
        color[u] = vcolor::black;
    };


    // DFS(G)を呼び出し、各頂点uに対して終了時刻u.fを計算する
    array_t tlst = tsort(G);

    // G^Tを計算する
    for (auto&& es : G) {
        for (auto&& e : es) {
            GT[e.dst].emplace_back(e.dst, e.src);
        }
    }

    // DFS(G^T)を呼び出す
    index_t k = 0;
    for (auto&& u : tlst) { // 成分グラフの頂点をトポロジカルソートされた順序で訪問する
        if (color[u] == vcolor::white) {
            dfs_visit(u, k++);
        }
    }

    // 分離された強連結成分を出力
    return components;
}



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



//****************************************
// エントリポイント
//****************************************

int main()
{
    using namespace std;
    using namespace graph;
    
    graph_t G;
    int s, t, m, n, q;
    cin >> n >> m;
    G.resize(n);
    for (index_t i = 0; i < m; i++) {
        cin >> s >> t;
        G[s].emplace_back(s, t);
        G[t].emplace_back(t, s);
    }

    auto _scc = scc(G);

    cin >> q;

    for (int i = 0; i < q; i++) {
        cin >> s >> t;
        cout << ((_scc[s] == _scc[t]) ? "yes" : "no") << endl;
    }

    return 0;
}

