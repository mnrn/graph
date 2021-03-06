/**
 * @brief グラフアルゴリズムに関する物置
 *
 * @note  グラフG = (V, E)には標準的な表現方法が2つある.隣接リストの集合による表現と隣接行列による表現である
 *        どちらを使っても有向グラフと無向グラフの両方を表現できる.グラフが疎(sparse)で|E|が|V|^2よりも
 *        ずっと小さいときには、隣接リストの表現はグラフをコンパクトに表現できる.これが表現を選択する際の基準になる
 *        グラフが密(dense)で|E|が|V|^2にほぼ等しいときや、2つの指定された頂点間に辺があるか否かを
 *        高速に判断する必要があるときには、隣接行列表現が好ましい場合もある
 *
 *        グラフG = (V, E)の隣接リスト表現(adjacency-list representation)は、Vの各頂点に対して1個、
 *        全部で|V|個のリストの配列Adjから構成される.各u ∈ Vに対して、隣接リストAdj[u]は辺(u, v) ∈ Eが存在する
 *        すべての頂点から構成される.(代わりに、隣接リストがこれらの頂点へのポインタを含むこともある.)
 *
 *        Gが有向グラフのとき、有向辺(u, v)はvが隣接リストAdj[u]の要素として現れることで、表現されるので、
 *        隣接リストの長さの合計は|E|である. Gが無向グラフのとき、(u, v)が無向辺ならば、uはvの隣接リストの要素であると
 *        同時にvはuの隣接リストの要素だから、隣接リストの長さの合計は2|E|である.有向グラフと無向グラフのどちらに対しても、
 *        隣接リスト表現はΘ(V + E)の記憶量しか必要としない望ましい性質を持っている
 *
 *        重み付きグラフ(weighted graph)では、重み関数(weight function) w:E -> R(Rは実数体)によって、各辺に対して
 *        その重みを定義する. 重み付きグラフを表現できるように隣接リストを変更するのは簡単である
 *        たとえば、G = (V, E)を重み関数wを持つ重み付きグラフとすると、辺(u, v) ∈ Eの重みw(u, v)をuの隣接リストに
 *        頂点vとともに格納するだけでよい. 多くのグラフの変形に提要できるように修正できるという点で、
 *        隣接リスト表現は非常に柔軟性に富んでいる
 *
 *        隣接リストが潜在的に持つ欠点は、与えられた辺(u, v)がグラフに属するか否かを決定するには、隣接リストAdj[i]の中から
 *        vを探索するより速い方法がないことである.この欠点は隣接行列表現を用いると救済できるが、
 *        それには漸近的により多くの記憶量が必要になる
 *
 *        グラフG = (V, E)の隣接行列表現(adjacency-matrix representation)では、ある方法で、番号0, 1,...,|V|-1が
 *        振られていると仮定する. グラフGの隣接行列表現は|V|x|V|型行列A = (aij)であり、
 *        aij = { 1 (i, j) ∈ Eのとき
 *                0 それ以外のとき }
 *        を満たす
 *        グラフの隣接行列はグラフの変数に関係なくΘ(V^2)の記憶量が必要である
 *
 *        隣接行列は主対角線に関して対称である. 無向グラフでは(u, v)と(v, u)は同じ辺を表すから、向こうグラフの隣接行列は
 *        その転置行列と等しい. すなわちA = Atである. いくつかの応用では、隣接行列の主対角線およびそれより上の要素だけを
 *        蓄え、グラフを格納するのに必要な記憶量をほぼ半分に減らすことができる
 *
 *        グラフの隣接リスト表現と同様、隣接行列表現でも重み付きグラフを表現できる. たとえば、 G = (V, E)を辺重み関数wを
 *        持つ重み付きグラフとするとき、辺(u, v) ∈ Eの重みw(u, v)を隣接行列のu行v列要素とすればよい
 *        辺が存在しないときは、NILを対応する行列要素の値としてもよいが、0や∞などの値を使ったほうが便利な場合も多い
 *
 *        隣接リスト表現の記憶効率は漸近的には隣接行列表現より悪くなることはないが、グラフが小さいときには
 *        隣接行列の単純さは捨てがたい. グラフが重み付きではない場合には1行列要素当たり1ビットで表現できるので、
 *        隣接行列表現はさらに有利になる
 *
 * @date  2016/02/08 ~ 2016/03/28
 */



//********************************************************************************
// インクルードガード
//********************************************************************************

 #ifndef GRAPH_HPP
 #define GRAPH_HPP



//********************************************************************************
// 必要なヘッダファイルのインクルード
//********************************************************************************

#include <vector>



//********************************************************************************
// オブジェクト形式マクロの定義
//********************************************************************************

#define GRAPH_BEGIN namespace graph {
#define GRAPH_END   }



//********************************************************************************
// 名前空間の始端
//********************************************************************************

GRAPH_BEGIN



//********************************************************************************
// 型シノニム
//********************************************************************************

using weight_t   = std::int32_t;  /**< 辺(u, v)への重みwを表す型 */
using index_t    = std::int32_t;  /**< 頂点vの添字を表す型       */
using capacity_t = weight_t;      /**< 辺(u, v)の容量を表す型    */
using bool_t     = std::int32_t;  /**< ブール値は整数型で、{ false 0 true それ以外 }とする */



//********************************************************************************
// 構造体の定義
//********************************************************************************

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
    vertex() noexcept : d(0), pi(0), color(vcolor::white)/*, f(0)*/ {}
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
    edge(index_t src, index_t dst) noexcept             : src(src), dst(dst), w(1) {}
    edge(index_t src, index_t dst, weight_t w) noexcept : src(src), dst(dst), w(w) {}
};



namespace limits {
    enum {  // scopedではあるが強く型付けされたenum(strongly-typed enum)ではない
        inf = std::numeric_limits<weight_t>::max() / 3,  /**< @brief 辺が存在しない場合に使用される値     */
        nil = std::numeric_limits<index_t>::min() / 3,   /**< @brief 先行点が存在しない場合に使用される値 */
    };
}



//********************************************************************************
// 型シノニムその2
//********************************************************************************

using edges_t    = std::vector<edge>;      /**< グラフG=(V, E)の辺集合E   */
using vertices_t = std::vector<vertex>;    /**< グラフG=(V, E)の頂点集合V */
using array_t    = std::vector<weight_t>;  /**< 重みwの配列  */
using indices_t  = std::vector<index_t>;   /**< 頂点の添字配列 */
using stamps_t   = std::vector<bool_t>;    /**< ブーリアンの集合 */
using matrix_t   = std::vector<array_t>;   /**< グラフGの隣接行列表現(および表行列表現) */
using graph_t    = std::vector<edges_t>;   /**< グラフGの隣接リスト表現(こちらを主に使用する) */



//********************************************************************************
// 名前空間の終端
//********************************************************************************

GRAPH_END



#endif  // end of GRAPH_HPP
