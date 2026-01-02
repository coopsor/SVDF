from __future__ import print_function

import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

import numpy as np
import tensorflow as tf

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(e)

from tensorflow.keras.layers import Flatten, Conv2D, MaxPooling2D
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Dropout, Input, Dense
from tensorflow.keras.layers import BatchNormalization


# from tensorflow.keras.utils.vis_utils import plot_model

def scaler_values(grouped_dataset):
    total_group = [list(j.data.values()) for i in grouped_dataset for j in i]
    X_max = np.max(total_group, axis=0)[:-1]
    X_min = np.min(total_group, axis=0)[:-1]
    scaler_dataset = []
    for group_data in grouped_dataset:
        intra_data = np.array([list(i.data.values()) for i in group_data])
        sublist = intra_data[:, :intra_data.shape[1] - 1]
        den = X_max - X_min
        den[den == 0] = 1
        scaler_sublist = (sublist - X_min) / den
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
               )(inputs)
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
    # plot_model(model, to_file='CNN_model.png', show_shapes=True)
    model.load_weights('model/model_class_test_015.h5')
    y_pred = np.argmax(model.predict(test_x), axis=1)

    filtered_clusters = []
    for pred_type, cluster_signature in zip(y_pred, grouped_dataset):
        if pred_type != 1:
            filtered_clusters.append(cluster_signature)

    return filtered_clusters
