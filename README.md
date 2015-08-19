indri-pruning-toolkit
====================

This toolkit implements various static index pruning algorithms with the Indri API.  

Currently, the following methods are supported:
- Term-based pruning (Carmel et al., 2001)
- Uniform pruning (Carmel et al., 2001; Chen and Lee, 2013)
- Document-centric pruning (Büttcher and Clarke, 2006)
- Probability ranking principle (Blanco and Barreiro, 2010)
- Two-sample two proportion test (Thota and Carterette, 2011)
- Popularity-based pruning (Ntoulas and Cho, 2007)
- Divergence based methods (Chen et al., 2015)

If you use this package in your research work, please cite the following paper:

> R.-C. Chen, C.-J. Lee, and W. B. Croft. On divergence measures and static index pruning. 
> In *Proceedings of ICTIR '15*. To appear.

### Installation ###

This package depends on:

1. An update-to-date version of g++ that supports C++0x, e.g., can do `-std=c++0x`;
2. Boost library >= 1.46;
3. A patched [Indri] library 5.x.

Better check if these requirements are satisfied prior to installation.  For
(3), you need to manually apply the patch at `share/DocListMemoryBuilder.patch`
to the Indri source.  

Suppose that the Indri source code (`indri-5.9` for example) and this repo
(`indri-pruning-toolkit`) are both placed under `$HOME/src`.  Do the following:

    cd $HOME/src/indri-5.9
    patch -p1 < $HOME/src/indri-pruning-toolkit/share/DocListMemoryBuilder.patch

Then compile the Indri code and get it installed.  Most likely you'll want to
put it under `$HOME`:

    ./configure --prefix=$HOME
    make -j8 && make install

If the first line doesn't work out and `configure` is really there, try `chmod
a+x configure`.  Once this comes through, go back to the toolkit repo and do:

    ./configure --prefix=$HOME --with-indri=$HOME
    make && make install


[Indri]: http://www.lemurproject.org/indri.php

### References ###

R. Blanco and A. Barreiro. Probabilistic static pruning of inverted files. *ACM
Trans. Inf. Syst.*, 28(1), Jan. 2010.

S. Büttcher and C. L. A. Clarke. A document-centric approach to static index
pruning in text retrieval systems. In *Proceedings of CIKM ’06*, pages 182–189.
ACM, 2006.

R.-C. Chen and C.-J. Lee. An information-theoretic account of static index
pruning.  In *Proceedings of SIGIR ’13*, pages 163–172. ACM, 2013.

R.-C. Chen, C.-J. Lee, and W. B. Croft. On divergence measures and static index
pruning. In *Proceedings of ICTIR '15*. To appear.

D. Carmel, D. Cohen, R. Fagin, E. Farchi, M. Herscovici, Y. S. Maarek, and A.
Soffer. Static index pruning for information retrieval systems. In *Proceedings
of SIGIR ’01*, pages 43–50. ACM, 2001.

A. Ntoulas and J. Cho. Pruning policies for two-tiered inverted index with
correctness guarantee. In *Proceedings of SIGIR ’07*, pages 191–198. ACM, 2007.

S. Thota and B. Carterette. Within-document term-based index pruning with
statistical hypothesis testing. In *Proceedings of ECIR ’11*, pages 543–554.
Springer Berlin Heidelberg, 2011.


