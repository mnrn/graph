/**
 * @brief  全点対最短路問題におけるFloyd-Warshallアルゴリズムを扱う
 * @note   関連URL: http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=GRL_1_C
 * @date   2016/02/21 ~ 2016/03/12
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
// 構造体の定義
//****************************************

using weight_t   = std::int32_t;  /**< 辺(u, v)への重みwを表す型 */
using index_t    = std::int32_t;  /**< 頂点vの添字を表す型       */
using capacity_t = weight_t;      /**< 辺(u, v)の容量を表す型    */



namespace limits {
    enum {  // scopedではあるが強く型付けされたenum(strongly-typed enum)ではない
        inf = std::numeric_limits<weight_t>::max(),  /**< @brief 辺が存在しない場合に使用される値     */
        nil = std::numeric_limits<index_t>::min(),   /**< @brief 先行点が存在しない場合に使用される値 */
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

/**
 * @brief  動的計画アルゴリズムであるFloyd-Warshallアルゴリズム
 *
 * @note   Floyd-Warshallアルゴリズムでは最短路の中間頂点を考える. ここで、単純道p = <v0,v1,...,vl>のv0とvl以外の頂点、
 *         すなわち、集合{ v1,v2,...,v(l-1) }の任意の頂点が中間(intermediate)頂点である
 *
 *         Floyd-Warshallアルゴリズムは以下の観察に基づいている. Gの頂点集合がV = { 0,1,...,n-1 }であるという仮定の下で、
 *         あるkに対して頂点の部分集合{ 0,1,...k }を考える. 任意の頂点対i,j ∈ Vに対して、iからjへの道で、そのすべての中間頂点が
 *         集合 { 0,1,...,k }に属するものをすべて考え、その中で最小の重みを持つ道をpとする. (道pは単純である)
 *
 *         Floyd-Warshallアルゴリズムは、道pと、中間頂点がすべて集合{ 0,1,...,k-1 }に属する(という条件の下での)iからjへの最短路との関係を解明する
 *         この関係はkが道pの中間頂点であるかどうかに依存する
 *
 *         * kが道pの中間頂点でなければ、pの中間頂点はすべて集合{ 0,1,...,k-1 }に属している. したがって、中間頂点がすべて集合{ 0,1,..,k-1 }に属する
 *           頂点iからjへの最短路は、中間頂点がすべて集合{ 0,1,...,k }に属するiからjへの最短路である
 *
 *         * kが道pの中間頂点ならば、pをi~(p1)~>k~(p2)~>jと分解する. p1は中間頂点がすべて集合{ 0,1,...,k }に属するiからjへの最短路である
 *           実際には、もう少し強い言明が成立する. 頂点kは道p1上の中間頂点ではないから、p1の中間頂点はすべて集合{ 0,1,...,k-1 }に属する.
 *           したがって、 p1は中間頂点がすべて集合{ 0,1,...,k-1 }に属するiからkへの最短路である
 *           同様に、p2は中間頂点がすべて集合{ 0,1,...,k-1 }に属する頂点kから頂点jへの最短路である
 *
 *         上の観察に基づいて、最短路の推定値を再帰的に定式化する. 中間頂点がすべて集合{ 0,...,k }に属する頂点iからjへの最短路の重みをdij^(k+1)とする
 *         k = -1のとき、-1より大きい番号を持つ中間頂点を通らない頂点iから頂点jへの道は中間頂点を1つも持たない. このような道は高々1本の辺しかもたないから、
 *         dij^(0) = wijである. よって、dij^(k+1)の再帰的な定義
 *           dij^(k+1) = { wij                            k = -1のとき
 *                       min(dij^(k), dik^(k) + djk^(k))  k >= 0のとき }
 *         を得る. 任意の道において、すべての中間頂点は集合{ 0,..,n-1 }に属するから、行列D^(n-1) = (dij^(n-1))が最終的な解を与える
 *         すなわち、すべてのi,j ∈ Vに対して、dij^(d) = δ(i, j)が成立する. 実際には、この場合を頂点0,1,...,k-1のみを使う場合に帰着させる
 *
 *         i,j,k = 0,1,...,n-1に対してdij^(k)を計算するからFloyd-WarshallアルゴリズムではΘ(n^3)の領域が必要であるが、
 *         前回の履歴のみを必要とするので、すべての肩添字を落としただけの手続きも正しく、Θ(n^2)の領域で十分である
 *
 * @note   3重のfor分における操作はΟ(1)時間で実行できるから、このアルゴリズムの実行時間はΘ(n^3)である
 *         このコードはタイトであって、複雑なデータ構造を含まず、Θ-記法に隠された定数は小さい
 *         したがって、Floyd-Warshallアルゴリズムは結構大きな入力グラフに対しても非常に実用的である
 *
 */
matrix_t floyd_warshall(const matrix_t& W)
{
    index_t n = W.size();  // n = W.rows
    matrix_t D(n, array_t(n, limits::inf));
    for (index_t i = 0; i < n; i++) { std::copy(W[i].begin(), W[i].end(), D[i].begin()); D[i][i] = 0; }

    
    for (index_t k = 0; k < n; k++) {
        for (index_t i = 0; i < n; i++) {
            for (index_t j = 0; j < n; j++) {
                // 任意の実数a != ∞に対して、a + ∞ = ∞ + a = ∞(また、a + (-∞) = (-∞) + a = (-∞))を仮定している
                if (D[i][k] != limits::inf && D[k][j] != limits::inf) {
                    D[i][j] = std::min(D[i][j], D[i][k] + D[k][j]);
                }
            }
        }
    }
    return D;
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
    
    int n, e, u, v, c;
    cin >> n >> e;
    matrix_t W(n, array_t(n, limits::inf));
    for (int i = 0; i < n; i++) { W[i][i] = 0; }

    for (int i = 0; i < e; i++) {
        cin >> u >> v >> c;
        W[u][v] = c;
    }

    matrix_t D = floyd_warshall(W);

    bool neg = false;
    for (int i = 0; i < n; i++) { if (D[i][i] < 0) { neg = true; } }

    if (neg) {
        cout << "NEGATIVE CYCLE" << endl;
    }
    else {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j) cout << " ";
                if (D[i][j] == limits::inf) cout << "INF";
                else                        cout << D[i][j];
            }
            cout << endl;
        }
    }
    
    return 0;
}

