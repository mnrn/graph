/**
 * @brief  最大フロー問題におけるFord-Fulkerson法を扱う
 * @note   関連URL: http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_6_A
 * @date   2016/02/24 ~ 2017/01/03
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>



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
// 型シノニム
//****************************************

using weight_t   = std::int32_t;  /**< 辺(u, v)への重みwを表す型 */
using index_t    = std::int32_t;  /**< 頂点vの添字を表す型       */
using capacity_t = weight_t;      /**< 辺(u, v)の容量を表す型    */
using bool_t     = std::int32_t;  /**< ブール値は整数型で、{ false 0 true それ以外 }とする */



//****************************************
// 構造体の定義
//****************************************

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
        bool_t  visited;  /**< 発見済みか?     */
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



namespace limits {
    enum {  // scopedではあるが強く型付けされたenum(strongly-typed enum)ではない
        inf = std::numeric_limits<weight_t>::max() / 3,  /**< @brief 辺が存在しない場合に使用される値     */
        nil = std::numeric_limits<index_t>::min() / 3,   /**< @brief 先行点が存在しない場合に使用される値 */
    };
}



//****************************************
// 型シノニム
//****************************************

using edges_t    = std::vector<edge>;      /**< グラフG=(V, E)の辺集合E   */
using vertices_t = std::vector<vertex>;    /**< グラフG=(V, E)の頂点集合V */
using array_t    = std::vector<weight_t>;  /**< 重みwの配列  */
using indices_t  = std::vector<index_t>;   /**< 頂点の添字配列 */
using stamps_t   = std::vector<bool_t>;    /**< ブーリアンの集合 */
using matrix_t   = std::vector<array_t>;   /**< グラフGの隣接行列表現(および表行列表現) */
using graph_t    = std::vector<edges_t>;   /**< グラフGの隣接リスト表現(こちらを主に使用する) */



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



//****************************************
// エントリポイント
//****************************************

int main(void)
{
    using namespace std;
    int V, E;

    cin >> V >> E;
    
    graph::graph_t G(V);
    for (int i = 0; i < E; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        G[a].emplace_back(a, b, c);
    }

    graph::ford_fulkerson ff(G);
    
    auto max_flow = ff.compute(0, V-1);

    cout << max_flow << endl;

    return 0;
}

