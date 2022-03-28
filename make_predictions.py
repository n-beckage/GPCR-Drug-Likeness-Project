import numpy as np
import sys
from multiprocessing import Pool
import os
import time
os.environ["CUDA_VISIBLE_DEVICES"]="0"    
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
import tensorflow as tf
import subprocess as sp
from tensorflow import keras
import subprocess as sp

# don't mess with this
def crossentropy(y_true,y_pred):
    return tf.reduce_mean(-1.*y_true[:,0]*tf.math.log(tf.clip_by_value(y_pred[:,0],1e-10,1.)) - (y_true[:,1])*tf.math.log(tf.clip_by_value(y_pred[:,1],1e-10,1.)))

# pass the whole array of fingerprints, and predictions will be made for each
def make_predictions(input_fpt): 
    activity_model=keras.models.load_model("model_cicular4.tf",custom_objects={'crossentropy':crossentropy}) # may need to adjust directory based on what the current wd is. An absolute path may be best
    activity_model.compile(
            loss=keras.losses.MeanSquaredError(),# keras.losses.MeanSquaredError(),
            optimizer=keras.optimizers.Adam(learning_rate=1),
            metrics=["MeanSquaredError"],
        )
    bsize=min([input_fpt.shape[0],1000]) # switch such that batch size does not exceed 1000
    if input_fpt.shape[0] !=0:
        pred_act=activity_model.predict(input_fpt,batch_size=bsize)
        return pred_act
    else:
        return 0


# Load model up here(first)

# Hint: ctrl + / will comment out highlighted text

fpts=np.load("fingerprints.npy",allow_pickle=True)
x=make_predictions(fpts)

# Last time I ran this script I got this error:
# ValueError: Failed to convert a NumPy array to a Tensor (Unsupported object type numpy.ndarray).