import keras
import numpy as np
import tables
from keras.models import Model
from keras.layers import Input, Dense, BatchNormalization

# 27 features
features = ['fj_jetNTracks',
            'fj_nSV',
            'fj_tau0_trackEtaRel_0',
            'fj_tau0_trackEtaRel_1',
            'fj_tau0_trackEtaRel_2',
            'fj_tau1_trackEtaRel_0',
            'fj_tau1_trackEtaRel_1',
            'fj_tau1_trackEtaRel_2',
            'fj_tau_flightDistance2dSig_0',
            'fj_tau_flightDistance2dSig_1',
            'fj_tau_vertexDeltaR_0',
            'fj_tau_vertexEnergyRatio_0',
            'fj_tau_vertexEnergyRatio_1',
            'fj_tau_vertexMass_0',
            'fj_tau_vertexMass_1',
            'fj_trackSip2dSigAboveBottom_0',
            'fj_trackSip2dSigAboveBottom_1',
            'fj_trackSip2dSigAboveCharm_0',
            'fj_trackSipdSig_0',
            'fj_trackSipdSig_0_0',
            'fj_trackSipdSig_0_1',
            'fj_trackSipdSig_1',
            'fj_trackSipdSig_1_0',
            'fj_trackSipdSig_1_1',
            'fj_trackSipdSig_2',
            'fj_trackSipdSig_3',
            'fj_z_ratio']

# 2 labels: QCD or Hbb
labels = ['fj_isQCD*sample_isQCD',
          'fj_isH*fj_isBB']

if __name__ == "__main__":
    """ This is executed when run from the command line """
    
    # load file
    h5file = tables.open_file('ntuple_merged_10.h5','r')
    njets = getattr(h5file.root,features[0]).shape[0]
    nfeatures = len(features)
    nlabels = len(labels)

    # allocate arrays
    feature_array = np.zeros((njets,nfeatures))
    label_array = np.zeros((njets,nlabels))

    # load feature arrays
    for (i, feat) in enumerate(features):
        feature_array[:,i] = getattr(h5file.root,feat)[:]

    # load labels arrays
    for (i, label) in enumerate(labels):
        prods = label.split('*')
        prod0 = prods[0]
        prod1 = prods[1]
        fact0 = getattr(h5file.root,prod0)[:]
        fact1 = getattr(h5file.root,prod1)[:]
        label_array[:,i] = np.multiply(fact0,fact1)

    # remove unlabeled samples
    print(label_array.shape)
    print(feature_array.shape)
    feature_array = feature_array[np.sum(label_array,axis=1)==1]
    label_array = label_array[np.sum(label_array,axis=1)==1]
    print(label_array.shape)
    print(feature_array.shape)

    # define dense keras model
    inputs = Input(shape=(nfeatures,), name = 'input')  
    x = BatchNormalization()(inputs)
    x = Dense(64, name = 'dense_1', activation='relu')(x)
    x = Dense(32, name = 'dense_2', activation='relu')(x)
    x = Dense(32, name = 'dense_3', activation='relu')(x)
    outputs = Dense(nlabels, name = 'output', activation='softmax')(x)
    keras_model = Model(inputs=inputs, outputs=outputs)
    keras_model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
    print(keras_model.summary())

    # fit keras model
    keras_model.fit(feature_array, label_array, batch_size=1024, 
                    epochs=100, validation_split=0.2, shuffle=False)
