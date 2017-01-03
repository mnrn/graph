# Graph Algorithm Library

## Overview
この業界に入る直前の期間に書いたグラフアルゴリズムプログラム集です。 
完全に学習用に書いたものなのでわりと独善的なコメントが目立ちます。 
ちょくちょく保守しています。


### Description

以下のアルゴリズムを扱っています。

* 基本グラフアルゴリズム(幅優先探索、深さ優先探索、トポロジカルソート、強連結成分)
* 最小全域木問題(KruskalとPrimのアルゴリズム)
* 単一始点最短路問題(Bellman-Ford法、Dijkstra法)
* 全点対最短路問題(Floyd-Warshall法)
* 最大フロー問題(Ford-Fulkerson法、Edmonds-Kerp法)


And open source with a [public repository][mnrn] on GitHub.

### Demo

[AOJ][AOJ]でアルゴリズムの動作確認をしました。 
各ディレクトリに存在するVerifyディレクトリに存在するsubmission.cppに確認した問題のURLが貼ってあります。


### Reference

以下のサイトをかなり参考にしました。
- [spaghetti source][spagetti source]
- [libalgo][libalgo]

基本的に[アルゴリズムイントロダクション][CLRS]の内容に沿ってプログラムを記述しています。


License
----

Public Domain


**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [mnrn]: <https://github.com/mnrn/graph>
   
   [spaghetti source]: <http://www.prefield.com/algorithm/>
   [libalgo]:<http://tubo28.me/algorithm/>

   [AOJ]: <http://judge.u-aizu.ac.jp/onlinejudge/>

   [CLRS]: <https://www.amazon.co.jp/%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0%E3%82%A4%E3%83%B3%E3%83%88%E3%83%AD%E3%83%80%E3%82%AF%E3%82%B7%E3%83%A7%E3%83%B3-%E7%AC%AC3%E7%89%88-%E7%B7%8F%E5%90%88%E7%89%88-%E4%B8%96%E7%95%8C%E6%A8%99%E6%BA%96MIT%E6%95%99%E7%A7%91%E6%9B%B8-%E3%82%B3%E3%83%AB%E3%83%A1%E3%83%B3/dp/476490408X>

