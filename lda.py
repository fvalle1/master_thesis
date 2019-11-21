import os
import sys
from typing import List, Any, Union

import numpy as np
from sklearn.decomposition import LatentDirichletAllocation
import pandas as pd
import tensorflow as tf
from sklearn.metrics import homogeneity_completeness_v_measure
from hsbmpy import get_file


def kullbach_liebler(theta_k, theta_l):
    # dey-visualizing paper
    return tf.subtract(tf.add(tf.math.multiply(theta_k, tf.math.log(tf.math.divide(theta_k, theta_l))), theta_l),
                       theta_k)


def distinctivness(theta_k):
    return tf.reduce_min(tf.sort(theta_k)[:, 1:], axis=1)


class lda(LatentDirichletAllocation):
    def __init__(self, learning_method='online',
                 max_doc_update_iter=5,
                 max_iter=5,
                 topic_word_prior=1,
                 doc_topic_prior=1,
                 random_state=42, **params):
        super().__init__(self, learning_method=learning_method, max_doc_update_iter=max_doc_update_iter,
                         topic_word_prior=topic_word_prior, doc_topic_prior=doc_topic_prior, max_iter=max_iter,
                         random_state=random_state, **params)
        if self.verbose > 1:
            print("model created")

    def __get_true_out(self, df, df_files):
        true_out = []
        for sample in df.columns.values:
            try:
                true_out.append(get_file(sample, df_files)['primary_site'])
            except:
                print(sys.exc_info()[1])
                true_out.append('')
        return true_out

    def full_analysis(self, directory, xl, tl=None, **kwargs) -> None:
        """

        :param df:
        :param directory:
        :param xl:
        :param tl:
        :param kwargs: argouments to LatentDirichletAllocation().fit_transform
        """
        scores = {}
        sigmas = []
        scores['lda'] = {
            'h': [],
            'c': [],
            'V': []
        }
        if tl is None:
            tl = xl
        df = pd.read_csv("%s/mainTable_all.csv" % directory, index_col=[0], header=[0]).dropna().astype(int)
        if self.verbose > 1:
            print(df.info())
        df_files = pd.read_csv("files.dat", index_col=[0])
        if self.verbose > 1:
            print(df_files.info())
        true_out=self.__get_true_out(df, df_files)

        if self.verbose > 1:
            print(df_files.info())
        totalobjcets = len(df.columns)
        print("lda")
        os.system('mkdir -p lda')
        for l, x in enumerate(xl):
            # lda
            ntopic = tl[l]
            # ntopic = x
            print("testing with %d clusters and %d topics" % (x, ntopic))
            self.n_components = ntopic
            print(self)
            topics = self.fit_transform(df.T.values, **kwargs)

            # save word distr
            print("saving word-distr")
            df_word_distr = pd.DataFrame(data=self.components_.T, index=df.index,
                                         columns=["Topic %d" % (t + 1) for t in np.arange(ntopic)])
            df_word_distr.to_csv("lda/lda_level_%d_word-dist.csv" % l, index=True, header=True)

            # save topic distr
            print("saving topic-distr")
            df_topic_distr = pd.DataFrame(data=topics, columns=["Topic %d" % (t + 1) for t in np.arange(ntopic)])
            df_topic_distr.insert(0, 'i_doc', np.arange(len(df.columns)))
            df_topic_distr.insert(1, 'doc', df.columns)
            df_topic_distr.to_csv("lda/lda_level_%d_topic-dist.csv" % l, index=False, header=True)
            del df_topic_distr

            # save topics
            data = tf.convert_to_tensor(df_word_distr.transpose())
            KL_tensor = tf.map_fn(fn=lambda k: tf.map_fn(fn=lambda l: kullbach_liebler(k, l), elems=data), elems=data)
            KL_tensor_min = tf.map_fn(distinctivness, tf.transpose(KL_tensor, perm=[2, 0, 1]))
            with tf.Session() as sess:
                out = sess.run(KL_tensor_min)
            df_D = pd.DataFrame(data=out, index=df_word_distr.index, columns=df_word_distr.columns)
            df_topics = pd.DataFrame(columns=df_D.columns, index=np.arange(len(df_D.index)))
            for topic in df_D.columns:
                df_topics[topic] = df_D[topic].sort_values(ascending=False).index
            df_topics.loc[:15, :].to_csv("%s/lda/lda_level_%d_topics.csv" % (directory, l), index=False)

            tf.reset_default_graph()

            # save clusters
            print("saving clusters")
            df_clusters = pd.DataFrame(index=np.arange(totalobjcets))
            # cluster = AgglomerativeClustering(n_clusters=x, affinity='euclidean', linkage='ward')
            # out = cluster.fit_predict(topics)
            out = np.argmax(topics, axis=1)

            for c in np.arange(out.max() + 1)[::-1]:
                c_objects = df.columns[np.argwhere(out == c)].values.T[0]
                df_clusters.insert(0, "Cluster %d" % (c + 1), np.concatenate(
                    (c_objects, [np.nan for _ in np.arange(totalobjcets - len(c_objects))])))
            df_clusters.dropna(axis=0, how='all', inplace=True)
            df_clusters.to_csv("lda/lda_level_%d_clusters.csv" % (l), index=False, header=True)

            # metrics
            print("saving metrics")
            score = (homogeneity_completeness_v_measure(true_out, out))
            scores['lda']['h'].append(score[0])
            scores['lda']['c'].append(score[1])
            scores['lda']['V'].append(score[2])

            # save dl
            sigmas.append(-self.score(X=df.values.T))

        pd.DataFrame(data=scores['lda']).to_csv("%s/lda.scores" % directory, header=True, index=False)


if __name__ == '__main__':
    model = lda(verbose=2, n_jobs=12)
    # model.full_analysis([[1,2,3],[4,5,6]])
