import pandas as pd
from scipy import stats
import numpy as np
from collections import Counter
import lightgbm as lgb
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFromModel
train_df = pd.read_csv("../data/train.csv")
test_df = pd.read_csv("../data/test.csv")

def choose_model_dat(model_type:str,train_data):
    '''
    model_type:str bact,fungi,metabolites,None
    '''
    if model_type=='bact':
        col_list=[i for i in train_data.columns if 'k__Bacteria' in i]
        # print(col_list)
        col_list=['Sample_ID']+col_list+['target']
        train_data=train_data[col_list]
    if model_type=='fungi':
        col_list=[i for i in train_data.columns if 'd__Fungi' in i]
        col_list=['Sample_ID']+col_list+['target']
        train_data=train_data[col_list]
    if model_type=='metabolites':
        col_list=[i for i in train_data.columns if ('d__Fungi' not in i) and ('k__Bacteria' not in i)]
        train_data=train_data[col_list]
    else:
        train_data=train_data
    return train_data


#---------------------------fungi---------------------------------------------------
train_df = pd.read_csv("../data/train.csv")
test_df = pd.read_csv("../data/test.csv")
train_df=choose_model_dat('fungi',train_df) 
train_df = train_df.drop(columns=['Sample_ID'])
test_df = test_df.drop(columns=['Sample_ID'])
print(train_df)

common_columns = train_df.columns.intersection(test_df.columns)
X_train = train_df[common_columns]
y_train = train_df['target']
X_test = test_df[common_columns]
y_test = test_df['M_stage']
missing_ratio = X_test.isnull().mean()
valid_columns = missing_ratio[missing_ratio <= 0.5].index
X_train = X_train[valid_columns]
X_test = X_test[valid_columns]
X_test = X_test.fillna(test_df.iloc[:,:-1].mean())
assert list(X_train.columns) == list(X_test.columns), "训练集和测试集的特征顺序不一致！"


col_lis_new = []
pi=0.1
for col in X_train.columns:
    group_0 = X_train[y_train == 'M0'][col]
    group_1 = X_train[y_train == 'M1'][col]
    p = stats.mannwhitneyu(group_0, group_1)
    if p.pvalue <= pi:
        col_lis_new.append(col)
print("Significant features:", col_lis_new)
X_train_selected = X_train[col_lis_new]
X_test_selected = X_test[col_lis_new]


X_data=X_train_selected
y_data=y_train

param_grid = {
    'num_leaves': [3,6,9,21],
    'learning_rate': [0.1,0.05,0.01,0.001],
    'n_estimators': [100,200],
    'max_depth': [3,4,5]
}

seed_range = range(1,51)
best_seed = None
best_params = None
best_score = -np.inf

results = []

for seed in seed_range:
    print(seed)
    X_train, X_test, y_train, y_test = train_test_split(
        X_data, y_data, test_size=0.2, random_state=seed, stratify=y_data
    )
    lasso = LogisticRegression(penalty='l1', solver='liblinear', random_state=42)
    lasso.fit(X_train, y_train)
    model = SelectFromModel(lasso, prefit=True)
    X_train_selected = model.transform(X_train)
    X_test_selected = model.transform(X_test)


    for num_leaves in param_grid['num_leaves']:
        for learning_rate in param_grid['learning_rate']:
            for n_estimators in param_grid['n_estimators']:
                for max_depth in param_grid['max_depth']:
                    
                    lgbm = lgb.LGBMClassifier(
                        objective='binary', 
                        metric='auc', 
                        random_state=seed,
                        device='gpu',
                        gpu_device_id=4,
                        num_leaves=num_leaves,
                        learning_rate=learning_rate,
                        n_estimators=n_estimators,
                        max_depth=max_depth
                    )
                    lgbm.fit(X_train_selected, y_train)

                    y_pred_prob = lgbm.predict_proba(X_test_selected)[:, 1]
                    score = roc_auc_score(y_test, y_pred_prob)
                    
                    results.append({
                        'seed': seed,
                        'num_leaves': num_leaves,
                        'learning_rate': learning_rate,
                        'n_estimators': n_estimators,
                        'max_depth': max_depth,
                        'score': score
                    })
                    
                    
                    if score > best_score:
                        best_score = score
                        best_seed = seed
                        best_params = {
                            'num_leaves': num_leaves,
                            'learning_rate': learning_rate,
                            'n_estimators': n_estimators,
                            'max_depth': max_depth
                        }


print(f"Best Seed: {best_seed}")
print(f"Best Parameters: {best_params}")
print(f"Best AUC Score: {best_score}")


with open(f"./result/M/best_result_{pi}_fungi.txt", "w") as f:
    f.write(f"Best Seed: {best_seed}\n")
    f.write(f"Best Parameters: {best_params}\n")
    f.write(f"Best AUC Score: {best_score:.4f}\n")


with open(f"./result/M/all_results.txt_{pi}_fungi", "w") as f:
    f.write("Seed\tNum_Leaves\tLearning_Rate\tN_Estimators\tMax_Depth\tScore\n")
    for res in results:
        f.write(f"{res['seed']}\t{res['num_leaves']}\t{res['learning_rate']}\t"
                f"{res['n_estimators']}\t{res['max_depth']}\t{res['score']:.4f}\n")     


sorted_results = sorted(results, key=lambda x: x['score'], reverse=True)
top_20_results = sorted_results[:30]       
with open(f"./result/M/top_20_results_{pi}_fungi.txt", "w") as f:
    f.write("Seed\tNum_Leaves\tLearning_Rate\tN_Estimators\tMax_Depth\tScore\n")
    for res in top_20_results:
        f.write(f"{res['seed']}\t{res['num_leaves']}\t{res['learning_rate']}\t"
                f"{res['n_estimators']}\t{res['max_depth']}\t{res['score']:.4f}\n")



#---------------------------bact---------------------------------------------------
train_df = pd.read_csv("../data/train.csv")
test_df = pd.read_csv("../data/test.csv")
train_df=choose_model_dat('bact',train_df) 
train_df = train_df.drop(columns=['Sample_ID'])
test_df = test_df.drop(columns=['Sample_ID'])
print(train_df)

common_columns = train_df.columns.intersection(test_df.columns)
X_train = train_df[common_columns]
y_train = train_df['target']
X_test = test_df[common_columns]
y_test = test_df['M_stage']
missing_ratio = X_test.isnull().mean()
valid_columns = missing_ratio[missing_ratio <= 0.5].index
X_train = X_train[valid_columns]
X_test = X_test[valid_columns]
X_test = X_test.fillna(test_df.iloc[:,:-1].mean())
assert list(X_train.columns) == list(X_test.columns), "训练集和测试集的特征顺序不一致！"


col_lis_new = []
pi=0.1
for col in X_train.columns:
    group_0 = X_train[y_train == 'M0'][col]
    group_1 = X_train[y_train == 'M1'][col]
    p = stats.mannwhitneyu(group_0, group_1)
    if p.pvalue <= pi:
        col_lis_new.append(col)
print("Significant features:", col_lis_new)
X_train_selected = X_train[col_lis_new]
X_test_selected = X_test[col_lis_new]


X_data=X_train_selected
y_data=y_train

param_grid = {
    'num_leaves': [3,6,9,21],
    'learning_rate': [0.1,0.05,0.01,0.001],
    'n_estimators': [100,200],
    'max_depth': [3,4,5]
}

seed_range = range(1,51)
best_seed = None
best_params = None
best_score = -np.inf

results = []

for seed in seed_range:
    print(seed)
    X_train, X_test, y_train, y_test = train_test_split(
        X_data, y_data, test_size=0.2, random_state=seed, stratify=y_data
    )
    lasso = LogisticRegression(penalty='l1', solver='liblinear', random_state=42)
    lasso.fit(X_train, y_train)
    model = SelectFromModel(lasso, prefit=True)
    X_train_selected = model.transform(X_train)
    X_test_selected = model.transform(X_test)


    for num_leaves in param_grid['num_leaves']:
        for learning_rate in param_grid['learning_rate']:
            for n_estimators in param_grid['n_estimators']:
                for max_depth in param_grid['max_depth']:
                    
                    lgbm = lgb.LGBMClassifier(
                        objective='binary', 
                        metric='auc', 
                        random_state=seed,
                        device='gpu',
                        gpu_device_id=4,
                        num_leaves=num_leaves,
                        learning_rate=learning_rate,
                        n_estimators=n_estimators,
                        max_depth=max_depth
                    )
                    lgbm.fit(X_train_selected, y_train)

                    y_pred_prob = lgbm.predict_proba(X_test_selected)[:, 1]
                    score = roc_auc_score(y_test, y_pred_prob)
                    
                    results.append({
                        'seed': seed,
                        'num_leaves': num_leaves,
                        'learning_rate': learning_rate,
                        'n_estimators': n_estimators,
                        'max_depth': max_depth,
                        'score': score
                    })
                    
                    
                    if score > best_score:
                        best_score = score
                        best_seed = seed
                        best_params = {
                            'num_leaves': num_leaves,
                            'learning_rate': learning_rate,
                            'n_estimators': n_estimators,
                            'max_depth': max_depth
                        }



print(f"Best Seed: {best_seed}")
print(f"Best Parameters: {best_params}")
print(f"Best AUC Score: {best_score}")


with open(f"./result/M/best_result_{pi}_bact.txt", "w") as f:
    f.write(f"Best Seed: {best_seed}\n")
    f.write(f"Best Parameters: {best_params}\n")
    f.write(f"Best AUC Score: {best_score:.4f}\n")


with open(f"./result/M/all_results.txt_{pi}_bact", "w") as f:
    f.write("Seed\tNum_Leaves\tLearning_Rate\tN_Estimators\tMax_Depth\tScore\n")
    for res in results:
        f.write(f"{res['seed']}\t{res['num_leaves']}\t{res['learning_rate']}\t"
                f"{res['n_estimators']}\t{res['max_depth']}\t{res['score']:.4f}\n")     


sorted_results = sorted(results, key=lambda x: x['score'], reverse=True)
top_20_results = sorted_results[:30]       
with open(f"./result/M/top_20_results_{pi}_bact.txt", "w") as f:
    f.write("Seed\tNum_Leaves\tLearning_Rate\tN_Estimators\tMax_Depth\tScore\n")
    for res in top_20_results:
        f.write(f"{res['seed']}\t{res['num_leaves']}\t{res['learning_rate']}\t"
                f"{res['n_estimators']}\t{res['max_depth']}\t{res['score']:.4f}\n")
        

#---------------------------metabolites---------------------------------------------------
train_df = pd.read_csv("../data/train.csv")
test_df = pd.read_csv("../data/test.csv")
train_df=choose_model_dat('metabolites',train_df) 
train_df = train_df.drop(columns=['Sample_ID'])
test_df = test_df.drop(columns=['Sample_ID'])
print(train_df)

common_columns = train_df.columns.intersection(test_df.columns)
X_train = train_df[common_columns]
y_train = train_df['target']
X_test = test_df[common_columns]
y_test = test_df['M_stage']
missing_ratio = X_test.isnull().mean()
valid_columns = missing_ratio[missing_ratio <= 0.5].index
X_train = X_train[valid_columns]
X_test = X_test[valid_columns]
X_test = X_test.fillna(test_df.iloc[:,:-1].mean())
assert list(X_train.columns) == list(X_test.columns), "训练集和测试集的特征顺序不一致！"


col_lis_new = []
pi=0.1
for col in X_train.columns:
    group_0 = X_train[y_train == 'M0'][col]
    group_1 = X_train[y_train == 'M1'][col]
    p = stats.mannwhitneyu(group_0, group_1)
    if p.pvalue <= pi:
        col_lis_new.append(col)
print("Significant features:", col_lis_new)
X_train_selected = X_train[col_lis_new]
X_test_selected = X_test[col_lis_new]


X_data=X_train_selected
y_data=y_train

param_grid = {
    'num_leaves': [3,6,9,21],
    'learning_rate': [0.1,0.05,0.01,0.001],
    'n_estimators': [100,200],
    'max_depth': [3,4,5]
}

seed_range = range(1,51)
best_seed = None
best_params = None
best_score = -np.inf

results = []

for seed in seed_range:
    print(seed)
    X_train, X_test, y_train, y_test = train_test_split(
        X_data, y_data, test_size=0.2, random_state=seed, stratify=y_data
    )
    lasso = LogisticRegression(penalty='l1', solver='liblinear', random_state=42)
    lasso.fit(X_train, y_train)
    model = SelectFromModel(lasso, prefit=True)
    X_train_selected = model.transform(X_train)
    X_test_selected = model.transform(X_test)


    for num_leaves in param_grid['num_leaves']:
        for learning_rate in param_grid['learning_rate']:
            for n_estimators in param_grid['n_estimators']:
                for max_depth in param_grid['max_depth']:
                    
                    lgbm = lgb.LGBMClassifier(
                        objective='binary', 
                        metric='auc', 
                        random_state=seed,
                        device='gpu',
                        gpu_device_id=4,
                        num_leaves=num_leaves,
                        learning_rate=learning_rate,
                        n_estimators=n_estimators,
                        max_depth=max_depth
                    )
                    lgbm.fit(X_train_selected, y_train)

                    y_pred_prob = lgbm.predict_proba(X_test_selected)[:, 1]
                    score = roc_auc_score(y_test, y_pred_prob)
                    
                    results.append({
                        'seed': seed,
                        'num_leaves': num_leaves,
                        'learning_rate': learning_rate,
                        'n_estimators': n_estimators,
                        'max_depth': max_depth,
                        'score': score
                    })
                    
                    
                    if score > best_score:
                        best_score = score
                        best_seed = seed
                        best_params = {
                            'num_leaves': num_leaves,
                            'learning_rate': learning_rate,
                            'n_estimators': n_estimators,
                            'max_depth': max_depth
                        }



print(f"Best Seed: {best_seed}")
print(f"Best Parameters: {best_params}")
print(f"Best AUC Score: {best_score}")


with open(f"./result/M/best_result_{pi}_metabolites.txt", "w") as f:
    f.write(f"Best Seed: {best_seed}\n")
    f.write(f"Best Parameters: {best_params}\n")
    f.write(f"Best AUC Score: {best_score:.4f}\n")


with open(f"./result/M/all_results.txt_{pi}_metabolites", "w") as f:
    f.write("Seed\tNum_Leaves\tLearning_Rate\tN_Estimators\tMax_Depth\tScore\n")
    for res in results:
        f.write(f"{res['seed']}\t{res['num_leaves']}\t{res['learning_rate']}\t"
                f"{res['n_estimators']}\t{res['max_depth']}\t{res['score']:.4f}\n")     


sorted_results = sorted(results, key=lambda x: x['score'], reverse=True)
top_20_results = sorted_results[:30]       
with open(f"./result/M/top_20_results_{pi}_metabolites.txt", "w") as f:
    f.write("Seed\tNum_Leaves\tLearning_Rate\tN_Estimators\tMax_Depth\tScore\n")
    for res in top_20_results:
        f.write(f"{res['seed']}\t{res['num_leaves']}\t{res['learning_rate']}\t"
                f"{res['n_estimators']}\t{res['max_depth']}\t{res['score']:.4f}\n")
        

#---------------------------all---------------------------------------------------
train_df = pd.read_csv("../data/train.csv")
test_df = pd.read_csv("../data/test.csv")
train_df=choose_model_dat('None',train_df) 
train_df = train_df.drop(columns=['Sample_ID'])
test_df = test_df.drop(columns=['Sample_ID'])
print(train_df)

common_columns = train_df.columns.intersection(test_df.columns)
X_train = train_df[common_columns]
y_train = train_df['target']
X_test = test_df[common_columns]
y_test = test_df['M_stage']
missing_ratio = X_test.isnull().mean()
valid_columns = missing_ratio[missing_ratio <= 0.5].index
X_train = X_train[valid_columns]
X_test = X_test[valid_columns]
X_test = X_test.fillna(test_df.iloc[:,:-1].mean())
assert list(X_train.columns) == list(X_test.columns), "训练集和测试集的特征顺序不一致！"


col_lis_new = []
pi=0.1
for col in X_train.columns:
    group_0 = X_train[y_train == 'M0'][col]
    group_1 = X_train[y_train == 'M1'][col]
    p = stats.mannwhitneyu(group_0, group_1)
    if p.pvalue <= pi:
        col_lis_new.append(col)
print("Significant features:", col_lis_new)
X_train_selected = X_train[col_lis_new]
X_test_selected = X_test[col_lis_new]


X_data=X_train_selected
y_data=y_train

param_grid = {
    'num_leaves': [3,6,9,21],
    'learning_rate': [0.1,0.05,0.01,0.001],
    'n_estimators': [100,200],
    'max_depth': [3,4,5]
}

seed_range = range(1,51)
best_seed = None
best_params = None
best_score = -np.inf

results = []

for seed in seed_range:
    print(seed)
    X_train, X_test, y_train, y_test = train_test_split(
        X_data, y_data, test_size=0.2, random_state=seed, stratify=y_data
    )
    lasso = LogisticRegression(penalty='l1', solver='liblinear', random_state=42)
    lasso.fit(X_train, y_train)
    model = SelectFromModel(lasso, prefit=True)
    X_train_selected = model.transform(X_train)
    X_test_selected = model.transform(X_test)


    for num_leaves in param_grid['num_leaves']:
        for learning_rate in param_grid['learning_rate']:
            for n_estimators in param_grid['n_estimators']:
                for max_depth in param_grid['max_depth']:
                    
                    lgbm = lgb.LGBMClassifier(
                        objective='binary', 
                        metric='auc', 
                        random_state=seed,
                        device='gpu',
                        gpu_device_id=4,
                        num_leaves=num_leaves,
                        learning_rate=learning_rate,
                        n_estimators=n_estimators,
                        max_depth=max_depth
                    )
                    lgbm.fit(X_train_selected, y_train)

                    y_pred_prob = lgbm.predict_proba(X_test_selected)[:, 1]
                    score = roc_auc_score(y_test, y_pred_prob)
                    
                    results.append({
                        'seed': seed,
                        'num_leaves': num_leaves,
                        'learning_rate': learning_rate,
                        'n_estimators': n_estimators,
                        'max_depth': max_depth,
                        'score': score
                    })
                    
                    if score > best_score:
                        best_score = score
                        best_seed = seed
                        best_params = {
                            'num_leaves': num_leaves,
                            'learning_rate': learning_rate,
                            'n_estimators': n_estimators,
                            'max_depth': max_depth
                        }


print(f"Best Seed: {best_seed}")
print(f"Best Parameters: {best_params}")
print(f"Best AUC Score: {best_score}")


with open(f"./result/M/best_result_{pi}_3model.txt", "w") as f:
    f.write(f"Best Seed: {best_seed}\n")
    f.write(f"Best Parameters: {best_params}\n")
    f.write(f"Best AUC Score: {best_score:.4f}\n")


with open(f"./result/M/all_results.txt_{pi}_3model", "w") as f:
    f.write("Seed\tNum_Leaves\tLearning_Rate\tN_Estimators\tMax_Depth\tScore\n")
    for res in results:
        f.write(f"{res['seed']}\t{res['num_leaves']}\t{res['learning_rate']}\t"
                f"{res['n_estimators']}\t{res['max_depth']}\t{res['score']:.4f}\n")     


sorted_results = sorted(results, key=lambda x: x['score'], reverse=True)
top_20_results = sorted_results[:30]       
with open(f"./result/M/top_20_results_{pi}_3model.txt", "w") as f:
    f.write("Seed\tNum_Leaves\tLearning_Rate\tN_Estimators\tMax_Depth\tScore\n")
    for res in top_20_results:
        f.write(f"{res['seed']}\t{res['num_leaves']}\t{res['learning_rate']}\t"
                f"{res['n_estimators']}\t{res['max_depth']}\t{res['score']:.4f}\n")