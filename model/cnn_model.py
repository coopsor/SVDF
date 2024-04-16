
import os
import pickle
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

import numpy as np
import tensorflow as tf
from tensorflow.python.keras.losses import categorical_crossentropy
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn import metrics
from tensorflow.python.keras.callbacks import LearningRateScheduler, EarlyStopping, ModelCheckpoint
from tensorflow.python.keras.layers import Flatten, Conv2D, MaxPooling2D
from tensorflow.python.keras.models import Model
from tensorflow.python.keras.layers import Dropout, BatchNormalization, Input, Dense
import tensorflow.python.keras.backend as K
from tensorflow.python.keras.backend import set_session

config = tf.compat.v1.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.3
set_session(tf.compat.v1.Session(config=config))

def to_one_hot(labels, dimension=3):
    results = np.zeros((len(labels), dimension))
    for i, label in enumerate(labels):
        results[i, int(label)] = 1.
    return results

def batch_generator(features, train_label, batch_size, max_length):
    # 计算需要多少批量
    num_batches = int(np.ceil(len(features) / batch_size))

    while True:
        for i in range(num_batches):
            # 获取当前批量的特征
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, len(features))
            batch_features = features[start_idx:end_idx]
            # 对每个特征进行填充，并将其合并为一个 numpy 数组
            batch_data = np.array([
                np.concatenate((matrix, np.zeros((max_length - matrix.shape[0], 25))), axis=0)
                for matrix in batch_features
            ])
            batch_data = np.expand_dims(batch_data, axis=-1)
            batch_label = train_label[start_idx:end_idx]
            batch_label = to_one_hot(batch_label)
            # batch_label = np.expand_dims(np.array(batch_label), axis=-1)
            yield batch_data, batch_label

train_X, train_Y, false_sam,false_sam_y,postive_sam,postive_sam_y = [], [], [],[],[],[]
with open('HG002/train_filter_aug.pkl', 'rb') as f:
    train_data = pickle.load(f)
    train_X, train_Y = train_data[0], train_data[1]

train_feature, val_feature, train_label, val_label = train_test_split(train_X, train_Y, test_size=0.3, random_state=43)
max_length = 100
verbose, epochs, batch_size = 1, 15, 256
n_timesteps, n_features, n_outputs = 100, 25, 1
checkpoint_path = "HG002/model/model_class_015.h5"

inputs = Input(shape=(max_length, n_features, 1))
x = Conv2D(filters=128, kernel_size=(3, n_features), padding='same', activation='relu',
           input_shape=(max_length, n_features, 1))(inputs)
x = MaxPooling2D(pool_size=(3, 1))(x)
x = Conv2D(filters=64, kernel_size=(3, 8), padding='same', activation='relu')(x)
x = MaxPooling2D(pool_size=(3, 1))(x)
x = Conv2D(filters=32, kernel_size=(3, 1), padding='same', activation='relu')(x)
x = MaxPooling2D(pool_size=(3, 1))(x)
# x = ECALayer()(x)
# x = BatchNormalization()(x)
x = Flatten()(x)
# x = RepeatVector(n_outputs)(x)
# x = Bidirectional(LSTM(256, activation='relu', return_sequences=True))(x)
# x = TimeDistributed(Dropout(0.5))(x)
# x = BatchNormalization()(x)
x = Dense(64, activation='relu')(x)
x = BatchNormalization()(x)
x = Dropout(0.5)(x)

x = Dense(32, activation='relu')(x)
x = BatchNormalization()(x)
x = Dropout(0.5)(x)
x = Dense(3, activation='softmax')(x)
model = Model(inputs, x)
model.compile(loss=categorical_crossentropy, optimizer='adam', metrics=['accuracy',])
model.summary()

def scheduler(epoch):
    if epoch % 2 == 0 and epoch != 0:
        lr = K.get_value(model.optimizer.lr)
        K.set_value(model.optimizer.lr, lr * 0.1)
        print("lr changed to {}".format(lr * 0.1))
    return K.get_value(model.optimizer.lr)

reduce_lr = LearningRateScheduler(scheduler)
early_stopping = EarlyStopping(patience=3, restore_best_weights=True)
train_gen = batch_generator(train_feature, train_label, batch_size, max_length)
valid_gen = batch_generator(val_feature, val_label, batch_size, max_length)

if True:
    model.fit(train_gen, steps_per_epoch=len(train_feature) // batch_size,
              epochs=epochs,
              batch_size=batch_size,
              verbose=verbose,
              validation_data=valid_gen,
              validation_steps=len(val_feature) // batch_size,
              callbacks=[reduce_lr,  ModelCheckpoint('HG002/model/model_class_test_{epoch:03d}.h5', save_freq='epoch')],
              )
else:
    model.load_weights(checkpoint_path)
test_x = np.array(
    [np.concatenate((matrix, np.zeros((max_length - matrix.shape[0], n_features))), axis=0) for matrix in
     val_feature])
test_x = np.reshape(test_x, (test_x.shape[0], test_x.shape[1], test_x.shape[2], 1))
y_pred = np.argmax(model.predict(test_x), axis=1)
print(metrics.classification_report(y_pred, val_label))
print('----------------------------------------------------------------------------------')

grouped_dataset = [i[0] for i in test_data]
label_list = [i[1] for i in test_data]
test_y = [i[0] for i in label_list]

test_x = np.array(
    [np.concatenate((matrix, np.zeros((max_length - matrix.shape[0], n_features))), axis=0) for matrix in
     grouped_dataset])
test_x = np.reshape(test_x, (test_x.shape[0], test_x.shape[1], test_x.shape[2], 1))

y_pred = np.argmax(model.predict(test_x), axis=1)

length_list = [len(i) for i in grouped_dataset]

# import matplotlib.pyplot as plt
# result = []
# for i, j, v in zip(y_pred, test_y, grouped_dataset):
#     if i != j:
#         result.append([j, len(v)])  # 真实标签，长度
# result = sorted(result, key=lambda x: x[1])
# result_array = np.array(result)
# first_data = result_array[:, 0]
# second_data = result_array[:, 1]
# categories = result_array[:, 0]
# counts = result_array[:, 1]
# # 统计第二维数据出现次数
# unique_counts, counts_occurrences = np.unique(counts, return_counts=True)
# # 绘制直方图
# plt.bar(unique_counts, counts_occurrences)
# plt.show()
#
print(metrics.classification_report(test_y, y_pred))

pred = []
for i, j in zip(length_list, y_pred):
    for k in range(i):
        pred.append(j)
true_label = [i[1] for i in label_list]
true_label = [j for i in true_label for j in i]
print(accuracy_score(true_label, pred))
print(metrics.classification_report(true_label, pred))
