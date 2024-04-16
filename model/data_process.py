import pickle
from collections import Counter
from multiprocessing import Pool
from random import sample

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.preprocessing import MinMaxScaler

filter_file = '../data/giab/HG002_SVs_Tier1_v0.6.bed'

def binary_search(ground_range, dataset, start, end):
    drop_index = []
    for index, row in enumerate(dataset[start:end]):
        row = list(row.data.values())
        # if row[-1] == 1:
        ground_list = ground_range[int(row[14])]
        target_start = row[6]
        target_end = row[7]
        flag = True
        for i in range(len(ground_list)):
            if ground_list[i][0] <= target_start and ground_list[i][1] >= target_end:
                flag = False
                break
        if flag:
            drop_index.append(index + start)

    return drop_index

def filter_data(dataset):
    filter_data = pd.read_csv(filter_file, sep='\t', header=None)
    chr_values = filter_data[0].values
    ground_range = []
    start_index = 0
    counter = Counter(chr_values)
    for element, count in counter.items():
        ground_range.append(filter_data.iloc[start_index:start_index + count, 1:3].values)
        start_index = start_index + count
    num_threads = 16
    chunk_size = len(dataset) // num_threads
    results = []
    analysis_pools = Pool(processes=int(num_threads))
    for i in range(num_threads):
        start = i * chunk_size
        end = start + chunk_size if i < num_threads - 1 else len(dataset)
        results.append(analysis_pools.apply_async(binary_search, (ground_range, dataset, start, end)))

    analysis_pools.close()
    analysis_pools.join()
    drop_index = [j for i in results for j in i.get()]
    dataset = np.delete(dataset, drop_index, axis=0)
    return dataset


# ont_50x = pd.read_csv('ont/intra_feature.txt', sep='\t', keep_default_na=False)
# ont_20x = pd.read_csv('ont_20x/intra_feature.txt', sep='\t', keep_default_na=False)
# ont_10x = pd.read_csv('ont_10x/intra_feature.txt', sep='\t', keep_default_na=False)
# ont_5x = pd.read_csv('ont_5x/intra_feature.txt', sep='\t', keep_default_na=False)
# ccs_30x = pd.read_csv('ccs/intra_feature.txt', sep='\t', keep_default_na=False)
# ccs_10x = pd.read_csv('ccs_10x/intra_feature.txt', sep='\t', keep_default_na=False)
# ccs_5x = pd.read_csv('ccs_5x/intra_feature.txt', sep='\t', keep_default_na=False)
# clr_65x = pd.read_csv('clr/intra_feature.txt', sep='\t', keep_default_na=False)
# clr_35x = pd.read_csv('clr_35x/intra_feature.txt', sep='\t', keep_default_na=False)
# clr_20x = pd.read_csv('clr_20x/intra_feature.txt', sep='\t', keep_default_na=False)
# clr_10x = pd.read_csv('clr_10x/intra_feature.txt', sep='\t', keep_default_na=False)
# clr_5x = pd.read_csv('clr_5x/intra_feature.txt', sep='\t', keep_default_na=False)

ont_signatures_50x = pickle.load(open('ont/signatures.pkl', 'rb'))[0]
ont_signatures_20x = pickle.load(open('ont_20x/signatures.pkl', 'rb'))[0]
ont_signatures_10x = pickle.load(open('ont_10x/signatures.pkl', 'rb'))[0]
ont_signatures_5x = pickle.load(open('ont_5x/signatures.pkl', 'rb'))[0]
ccs_signatures_30x = pickle.load(open('ccs/signatures.pkl', 'rb'))[0]
ccs_signatures_10x = pickle.load(open('ccs_10x/signatures.pkl', 'rb'))[0]
ccs_signatures_5x = pickle.load(open('ccs_5x/signatures.pkl', 'rb'))[0]
clr_signatures_65x = pickle.load(open('clr/signatures.pkl', 'rb'))[0]
clr_signatures_35x = pickle.load(open('clr_35x/signatures.pkl', 'rb'))[0]
clr_signatures_20x = pickle.load(open('clr_20x/signatures.pkl', 'rb'))
clr_signatures_10x = pickle.load(open('clr_10x/signatures.pkl', 'rb'))
clr_signatures_5x = pickle.load(open('clr_5x/signatures.pkl', 'rb'))[0]

def form_partitions(sv_signatures, max_distance):
    """Form partitions of signatures using mean distance."""
    sorted_signatures = sorted(sv_signatures, key=lambda evi: evi.get_key())
    partitions = []
    current_partition = []
    for signature in sorted_signatures:
        if len(current_partition) > 0 and current_partition[-1].downstream_distance_to(signature) > max_distance:
            partitions.append(current_partition[:])
            current_partition = []
        current_partition.append(signature)
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    return partitions

def span_position_distance(signature1, signature2, type, param):
    span1 = signature1.get_source()[2] - signature1.get_source()[1]
    span2 = signature2.get_source()[2] - signature2.get_source()[1]
    center1 = (signature1.get_source()[1] + signature1.get_source()[2]) // 2
    center2 = (signature2.get_source()[1] + signature2.get_source()[2]) // 2
    span_distance = abs(span1 - span2) / max(span1, span2)
    postion_distance = abs(center1 - center2)
    return postion_distance // (max(span1, span2) * param) + span_distance

def group_data(raw_data, temp):
    partitions = form_partitions(raw_data, 1000)
    clusters_final = []
    large_partitions = 0
    # Find clusters in each partition individually.
    cluster_len = [len(partition) for partition in partitions]
    mean_len = sum(cluster_len) / len(cluster_len)
    print(mean_len)
    for partition in partitions:
        if len(partition) == 1:
            if temp in [4, 5, 6]:
                clusters_final.append(partition)
            continue
        elif len(partition) > 100:
            partition_sample = sample(partition, 100)
            partition_sample.sort(key=lambda signature: (signature.contig, signature.start, signature.end))
            large_partitions += 1
        else:
            partition_sample = partition
        element_type = partition_sample[0].type

        param = (abs(len(partition_sample) - mean_len)) / max(len(partition_sample), mean_len) + mean_len
        distance_data = []
        for i in range(len(partition_sample) - 1):  # 局部深度反映了当前sv在1000bp以内的相似sv的个数
            for j in range(i + 1, len(partition_sample)):
                distance_data.append(
                    span_position_distance(partition_sample[i], partition_sample[j], element_type, param))
        Z = linkage(np.array(distance_data), method="average")
        cluster_max = 0.3
        cluster_indices = list(fcluster(Z, cluster_max, criterion='distance'))
        new_clusters = [[] for i in range(max(cluster_indices))]
        for signature_index, cluster_index in enumerate(cluster_indices):
            new_clusters[cluster_index - 1].append(partition_sample[signature_index])
        clusters_final.extend(new_clusters)

    return clusters_final

max_length = 100
scaler = MinMaxScaler()

def scaler_test_values(grouped_dataset, flag):
    total_group = [list(j.data.values()) for i in grouped_dataset for j in i]
    X_max = np.max(total_group, axis=0)[:-1]
    X_min = np.min(total_group, axis=0)[:-1]
    scaler_dataset = []
    for group_data in grouped_dataset:
        intra_data = np.array([list(i.data.values()) for i in group_data])
        sublist = intra_data[:, :intra_data.shape[1] - 1]
        label = intra_data[:, -1].tolist()
        for i, j in enumerate(label):
            if j != 1:
                label[i] = 0 if sublist[i][11] == 1 else 2
        # label = [0 if i != 1 else 1 for i in label]
        label_y = (max(label, key=label.count))
        scaler_sublist = (sublist - X_min) / (X_max - X_min)
        scaler_sublist = np.concatenate((scaler_sublist, np.tile(flag, (scaler_sublist.shape[0], 1))), axis=1)
        if len(sublist) <= max_length:
            scaler_dataset.append([scaler_sublist, (label_y, label)])
        elif len(sublist) > max_length:
            for i in range(0, len(sublist) // max_length):
                start = i * max_length
                end = (i + 1) * max_length if (i + 1) * max_length < len(sublist) else len(sublist)
                scaler_dataset.append([scaler_sublist[start:end], (label_y, label[start:end])])

    return scaler_dataset

def augument_data(train_X, train_Y):
    # 数据增强
    augmented_data = []
    augmented_data_label = []
    for x, y in zip(train_X, train_Y):
        if x.shape[0] > 20:
            continue
        for target_length in range(2, 20):
            for i in range(x.shape[0]):
                if x.shape[0] - target_length <= 0:
                    break
                starts = np.random.randint(low=0, high=x.shape[0] - target_length)
                ends = target_length + starts
                windows_slice_data = x[starts:ends]
                augmented_data.append(windows_slice_data)
                augmented_data_label.append(y)
        augmented_data.append(x)
        augmented_data_label.append(y)
    element_counts = Counter(augmented_data_label)
    for element, count in element_counts.items():
        print(f"数据增强后元素 {element} 出现了 {count} 次")
    return augmented_data, augmented_data_label

train_grouped_dataset, test_dataset, test_grouped_dataset = [], [], []
train_data, train_label, test_data = [], [], []
temp = 0

for sub_signature in [clr_signatures_65x, clr_signatures_35x, clr_signatures_20x, clr_signatures_10x, clr_signatures_5x,
                      ccs_signatures_30x, ccs_signatures_10x, ccs_signatures_5x, ont_signatures_50x, ont_signatures_20x,
                      ont_signatures_10x, ont_signatures_5x]:
    temp += 1
    sub_train_dataset, sub_train_grouped = [], []
    for i in sub_signature:
        if not hasattr(i, 'data'):
            continue
        if i.data['C'] < 11:
            sub_train_dataset.append(i)
    sub_train_dataset = filter_data(sub_train_dataset)
    with open('HG002/train_data_{}.pkl'.format(temp), 'wb') as f:
        pickle.dump(sub_train_dataset, f)

for temp in range(12):
    if temp in [0, 1, 2, 3]:
        flag = 0
    elif temp in [4, 5, 6]:
        flag = 1
    elif temp in [7, 8, 9, 10, 11]:
        flag = 2
    plat = [0, 0, 0]
    plat[flag] = 1
    with open('HG002/train_data_{}.pkl'.format(temp), 'rb') as f:
        sub_train_dataset = pickle.load(f)
    del_sig = [i for i in sub_train_dataset if i.type == 'DEL']
    ins_sig = [i for i in sub_train_dataset if i.type == 'INS']
    sub_train_grouped = group_data(del_sig, temp)
    sub_ins_grouped = group_data(ins_sig, temp)
    sub_train_grouped.extend(sub_ins_grouped)
    scaler_train_data = scaler_test_values(sub_train_grouped, plat)
    train_X = [i[0] for i in scaler_train_data]
    label_list = [i[1] for i in scaler_train_data]
    train_Y = [i[0] for i in label_list]
    element_counts = Counter(train_Y)
    for element, count in element_counts.items():
        print(f"元素 {element} 出现了 {count} 次")
    train_X, train_Y = augument_data(train_X, train_Y)
    train_data.extend(train_X)
    train_label.extend(train_Y)
    temp += 1

element_counts = Counter(train_label)
for element, count in element_counts.items():
    print(f"元素 {element} 出现了 {count} 次")
with open('HG002/train_filter_aug.pkl', 'wb') as f:
    pickle.dump([train_data, train_label], f)

# lengths = [[len(sublist), label] for sublist, label in zip(train_feature, train_y)]
# lengths.sort(reverse=True)
# import matplotlib.pyplot as plt
#
# values_0 = [item[0] for item in lengths if item[1] == 0]
# values_1 = [item[0] for item in lengths if item[1] == 1]
#
# hist_0, bins_0 = np.histogram(values_0, bins=np.unique(values_0))
# hist_1, bins_1 = np.histogram(values_1, bins=np.unique(values_1))
#
# plt.bar(bins_0[:-1], hist_0, width=0.5, align='center', color='blue', alpha=0.7, label='Label 0')
# for value, freq in zip(bins_0[:-1], hist_0):
#     plt.text(value, freq, f'({value},{freq})', ha='center', va='bottom')
# plt.show()
# plt.bar(bins_1[:-1], hist_1, width=0.5, align='center', color='green', alpha=0.7, label='Label 1')
#
# for value, freq in zip(bins_1[:-1], hist_1):
#     plt.text(value, freq, f'({value},{freq})', ha='center', va='bottom')
# plt.show()
# #

# data_see = []
# data_around = []
# i = 0
# for sublist, label in zip(train_feature, train_y):
#     if i - 2 > 0 and i + 2 < len(train_feature) and label == 1:
#         if len(sublist) == 1 or len(sublist) == 2:
#             data_see.append([sub_sublist[5] for sub_sublist in sublist])
#             data_around.append([sub_sublist for sub_sublist in train_feature[i-2:i+2]])
#     i = i+1


# train_X, train_Y = augument_data()
# train_feature, val_feature, train_label, val_label = train_test_split(train_X, train_Y, test_size=0.1, random_state=42)


# lengths = [[len(sublist), label] for sublist, label in zip(train_feature, train_y)]
# lengths.sort(reverse=True)
#
# import matplotlib.pyplot as plt
#
# values_0 = [item[0] for item in lengths if item[1] == 0]
# values_1 = [item[0] for item in lengths if item[1] == 1]
#
# hist_0, bins_0 = np.histogram(values_0, bins=np.unique(values_0))
# hist_1, bins_1 = np.histogram(values_1, bins=np.unique(values_1))
#
# plt.bar(bins_0[:-1], hist_0, width=0.5, align='center', color='blue', alpha=0.7, label='Label 0')
# for value, freq in zip(bins_0[:-1], hist_0):
#     plt.text(value, freq, f'({value},{freq})', ha='center', va='bottom')
# plt.show()
# plt.bar(bins_1[:-1], hist_1, width=0.5, align='center', color='green', alpha=0.7, label='Label 1')
#
# for value, freq in zip(bins_1[:-1], hist_1):
#     plt.text(value, freq, f'({value},{freq})', ha='center', va='bottom')
# plt.show()
