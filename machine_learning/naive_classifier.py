import numpy as np

class NaiveClassifier:
    
    #average max score of positive samples in training data
    _avg_max_score = 0

    #maximum max score of positive samples in training data
    _max_max = 0

    def fit(self, train_features, train_label):
        #get all indices with 1 in train_label
        indices = []
        for ind, label in enumerate(train_label):
            if label == 1:
                indices.append(ind)
        
        #get the max chaining score of all positive samples
        pos_samples = train_features.iloc[indices]
        pos_samples_max = pos_samples.max(axis="columns")
        
        #find average max
        sum = 0
        num_items = 0
        #find maximum max
        max_max = 0
        for item in pos_samples_max:
            sum += item
            num_items += 1
            if item > max_max:
                max_max = item
        
        self._avg_max_score = sum / num_items
        self._max_max = max_max
        return

    def predict(self, test_features):
        samples_max = test_features.max(axis="columns")
        predictions = []
        for item in samples_max:
            if item >= self._avg_max_score:
                predictions.append(1)
            else:
                predictions.append(0)
        return predictions
    
    '''
    Returns an array of shape (n_samples, n_classes). Each element of the array corresponds to
    probability that the corresponding sample belongs to the corresponding class. Classes
    are in the order [0, 1].
    '''
    def predict_proba(self, test_features):
        samples_max = test_features.max(axis="columns")
        predictions = []
        for item in samples_max:
            prob_1 = item / self._max_max
            if prob_1 > 1:
                prob_1 = 1
            prob_0 = 1 - prob_1
            predictions.append([prob_0, prob_1])
        return np.array(predictions)