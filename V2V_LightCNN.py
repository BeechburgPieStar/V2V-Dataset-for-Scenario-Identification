# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 14:50:49 2020

@author: Noah
"""

import scipy.io as scio
import numpy as np
# import pandas as pd
from keras.models import Model
from sklearn.metrics import confusion_matrix
from numpy import array
from keras.utils import to_categorical
from sklearn.model_selection import train_test_split
import keras.models as models
from keras.layers import *
from keras.callbacks import ModelCheckpoint
from keras import regularizers

KS = 3
def CNN():
	cnn_input = Input(shape=[52, 2])
	x = Conv1D(32, KS, padding='same')(cnn_input)
	x = Activation('relu')(x)
	x = SeparableConv1D(32, KS, padding='same')(x)
	x = Activation('relu')(x)

	x = AveragePooling1D(4)(x)

	x = SeparableConv1D(32, KS, padding='same')(x)
	x = Activation('relu')(x)
	x = SeparableConv1D(32, KS, padding='same')(x)
	x = Activation('relu')(x)

	x = GlobalAveragePooling1D()(x)

	x = Dense(5, activation="softmax", kernel_regularizer=regularizers.l2(0.01))(x)

	model = Model(inputs=cnn_input, outputs=x)
	model.compile(loss='categorical_crossentropy',optimizer='adam', metrics=['accuracy'])
	return model

model = CNN()
model.summary()

data_path="V2V_H.mat"
data = scio.loadmat(data_path)
x=data.get('H')
Sample_num = 10000
y1=np.zeros([Sample_num,1])
y2=np.ones([Sample_num,1])
y3=np.ones([Sample_num,1])*2
y4=np.ones([Sample_num,1])*3
y5=np.ones([Sample_num,1])*4
y=np.vstack((y1,y2,y3,y4,y5))
y = array(y)
y = to_categorical(y)
X_train, X_val, Y_train, Y_val = train_test_split(x, y, test_size = 0.3, random_state= 30)
model = CNN()
model.summary()
checkpoint = ModelCheckpoint("LCNN_GAP.hdf5", verbose=1, save_best_only=True)

hist=model.fit(
    X_train,
    Y_train,
    batch_size=100,
    epochs=100,
    verbose=1,
    validation_data=(X_val, Y_val),
    callbacks=[checkpoint]
    )
train_test_list = [hist.history['acc'],hist.history['val_acc'],hist.history['loss'],hist.history['val_loss']]
train_test_array=np.array(train_test_list).T
df = pd.DataFrame(train_test_array, columns=['Training Acc', 'Test Acc','Training Loss','Test Loss'])
df.to_excel("LCNN_GAP.xlsx", index=False)

import time
SNRs = range(15, 41, 1)
for snr in SNRs:
    data_path="V2V_H_SNR="+str(snr)+".mat"
    data = scio.loadmat(data_path)
    x=data.get('H')
    Sample_num = 10000
    y1=np.zeros([Sample_num,1])
    y2=np.ones([Sample_num,1])
    y3=np.ones([Sample_num,1])*2
    y4=np.ones([Sample_num,1])*3
    y5=np.ones([Sample_num,1])*4
    y=np.vstack((y1,y2,y3,y4,y5))
    y = array(y)
    y = to_categorical(y)
    # model.load_weights("LCNN_GAP.hdf5")
    t1 = time.time()
    [loss, acc] = model.evaluate(x, y, batch_size = 100, verbose=0)
    t2 = time.time()
    print(t2-t1) 