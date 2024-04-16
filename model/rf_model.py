import pickle

import pandas as pd
import sklearn.metrics as metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split

with open('pbsim_ont_v2/sig_inter_.pkl', 'rb') as temp:
    sig_inter = pickle.load(temp)
with open('pbsim_clr_v2/sig_inter_.pkl', 'rb') as temp:
    sig_inter_clr = pickle.load(temp)
df_indel = pd.read_csv('inter_data/train_inter_ont.txt', sep='\t',  keep_default_na=False)
df_indel = df_indel[df_indel['label'] != 'None']
df_indel_2 = pd.read_csv('inter_data/train_inter_clr.txt', sep='\t',  keep_default_na=False)
df_indel_2 = df_indel_2[df_indel_2['label'] != 'None']
df_indel_3 = pd.read_csv('inter_data/train_inter.txt', sep='\t',  keep_default_na=False)
#将最后一列中DEL,INS,None分别替换为0,1,2
for df in [df_indel, df_indel_2, df_indel_3]:
    df['label'] = df['label'].replace('DEL', 0)
    df['label'] = df['label'].replace('INS', 1)
    df['label'] = df['label'].replace('None', 5)

data = []
for sv in sig_inter:
    data.append(sv.data[0])
for sv in sig_inter_clr:
    data.append(sv.data[0])
df_ = pd.DataFrame(data)
df_ = pd.concat([pd.DataFrame(data), df_indel, df_indel_2, df_indel_3])
df_ = df_.drop_duplicates()
df_ = df_.dropna(axis=0, how='any')

df = df_[df_['C']<11]
X = df.iloc[:, 0:-1]
# X = StandardScaler().fit_transform(X.values)
Y =df['label']
print(Y.value_counts())
train_X, test_X, train_Y, test_Y = train_test_split(X.values, Y.values, test_size=0.3, stratify=Y,
                                                    random_state=42)
model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')

model.fit(train_X, train_Y)
with open('model_sim_indel_rf_test_.pkl', 'wb') as f:
    pickle.dump(model, f)
with open('model_sim_indel_rf_test_.pkl', 'rb') as f:
    model = pickle.load(f)
pred = model.predict(test_X)
print(pred)
# 计算f1得分
f1 = f1_score(test_Y, pred, average='weighted')
print(f'F1 score: {f1}')
print(metrics.classification_report(test_Y, pred))

importances = model.feature_importances_
# 打印每个特征的重要性
for i, importance in enumerate(importances):
    print(f"Feature {i}: {importance}")
