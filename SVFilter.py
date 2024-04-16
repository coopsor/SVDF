from __future__ import print_function

import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

import numpy as np
import tensorflow as tf
from tensorflow.python.keras.backend import set_session
config = tf.compat.v1.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.3
set_session(tf.compat.v1.Session(config=config))

from tensorflow.python.keras.layers import Flatten, Conv2D, MaxPooling2D
from tensorflow.python.keras.models import Model
from tensorflow.python.keras.layers import Dropout, Input, Dense
from tensorflow.python.keras.layers import BatchNormalization
from tensorflow.python.keras.utils.vis_utils import plot_model

def scaler_values(grouped_dataset):
    total_group = [list(j.data.values()) for i in grouped_dataset for j in i]
    X_max = np.max(total_group, axis=0)[:-1]
    X_min = np.min(total_group, axis=0)[:-1]
    scaler_dataset = []
    for group_data in grouped_dataset:
        intra_data = np.array([list(i.data.values()) for i in group_data])
        sublist = intra_data[:, :intra_data.shape[1] - 1]
        scaler_sublist = (sublist - X_min) / (X_max - X_min)
        scaler_dataset.append(scaler_sublist)
    return scaler_dataset

def filter_clusters(grouped_dataset):
    max_length, n_features, n_outputs = 100, 25, 1
    scaler_dataset = scaler_values(grouped_dataset)
    test_x = np.array(
        [np.concatenate((matrix, np.zeros((max_length - matrix.shape[0], n_features))), axis=0) for matrix in
         scaler_dataset])
    test_x = np.reshape(test_x, (test_x.shape[0], test_x.shape[1], test_x.shape[2], 1))

    inputs = Input(shape=(max_length, n_features, 1))
    x = Conv2D(filters=128, kernel_size=(3, n_features), padding='same', activation='relu',
               input_shape=(50, n_features, 1))(inputs)
    x = MaxPooling2D(pool_size=(3, 1))(x)
    x = Conv2D(filters=64, kernel_size=(3, 8), padding='same', activation='relu')(x)
    x = MaxPooling2D(pool_size=(3, 1))(x)
    x = Conv2D(filters=32, kernel_size=(3, 1), padding='same', activation='relu')(x)
    x = MaxPooling2D(pool_size=(3, 1))(x)
    x = Flatten()(x)
    x = Dense(64, activation='relu')(x)
    x = BatchNormalization()(x)
    x = Dropout(0.5)(x)
    x = Dense(32, activation='relu')(x)
    x = BatchNormalization()(x)
    x = Dropout(0.5)(x)
    x = Dense(3, activation='softmax')(x)
    model = Model(inputs, x)
    plot_model(model, to_file='CNN_model.png', show_shapes=True)
    model.load_weights('model/model_class_test_015.h5')
    y_pred = np.argmax(model.predict(test_x), axis=1)
    from numba import cuda
    cuda.select_device(0)
    cuda.close()
    filtered_clusters = []
    for pred_type, cluster_signature in zip(y_pred, grouped_dataset):
        if pred_type != 1:
            filtered_clusters.append(cluster_signature)
            
    return filtered_clusters
