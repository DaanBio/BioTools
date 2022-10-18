import io
import os
import scipy.misc
import numpy as np
import six
import time
import glob
from IPython.display import display

from six import BytesIO

import matplotlib
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

import tensorflow as tf
from object_detection.utils import ops as utils_ops
from object_detection.utils import label_map_util
from object_detection.utils import visualization_utils as vis_util
import streamlit as st

from tqdm import tqdm

# streamlit part
st.title("Colony finder")
st.markdown('This app counts CFUs and shows the location of the CFUs on a petridish. The app is meant to help students saving time counting colonies'
            ' or determining sizes of their colonies. The app uses a Region based Convolutional Neural Network based on efficientnet d1.')
uploaded_file = st.file_uploader("Upload Files", type=['png', 'jpeg','jpg'])

labelmap_path='Colonyfinder/label_map.pbtxt'
category_index = label_map_util.create_category_index_from_labelmap(labelmap_path, use_display_name=True)

with st.spinner('loading R-CNN model....'):
    tf.keras.backend.clear_session()
    model = tf.saved_model.load('Colonyfinder/saved_model')
st.success('Done!')

def run_inference_for_single_image(model, image):
    image = np.asarray(image)
    # The input needs to be a tensor, convert it using `tf.convert_to_tensor`.
    input_tensor = tf.convert_to_tensor(image)
    # The model expects a batch of images, so add an axis with `tf.newaxis`.
    input_tensor = input_tensor[tf.newaxis, ...]

    # Run inference
    model_fn = model.signatures['serving_default']
    output_dict = model_fn(input_tensor)

    # All outputs are batches tensors.
    # Convert to numpy arrays, and take index [0] to remove the batch dimension.
    # We're only interested in the first num_detections.
    num_detections = int(output_dict.pop('num_detections'))
    output_dict = {key: value[0, :num_detections].numpy()
                   for key, value in output_dict.items()}
    output_dict['num_detections'] = num_detections

    # detection_classes should be ints.
    output_dict['detection_classes'] = output_dict['detection_classes'].astype(np.int64)

    # Handle models with masks:
    if 'detection_masks' in output_dict:
        # Reframe the the bbox mask to the image size.
        detection_masks_reframed = utils_ops.reframe_box_masks_to_image_masks(
            output_dict['detection_masks'], output_dict['detection_boxes'],
            image.shape[0], image.shape[1])
        detection_masks_reframed = tf.cast(detection_masks_reframed > 0.5,
                                           tf.uint8)
        output_dict['detection_masks_reframed'] = detection_masks_reframed.numpy()

    return output_dict
if uploaded_file is not None:
    with st.spinner('Running inference on the image....'):

        img = Image.open(uploaded_file)
        image_np = np.array(img)
        output_dict = run_inference_for_single_image(model, image_np)
        vis_util.visualize_boxes_and_labels_on_image_array(
              image_np,
              output_dict['detection_boxes'],
              output_dict['detection_classes'],
              output_dict['detection_scores'],
              category_index,
              instance_masks=output_dict.get('detection_masks_reframed', None),
              use_normalized_coordinates=True,
              max_boxes_to_draw=100,
              line_thickness=12)
        check = (output_dict['detection_scores'] >= 0.5).sum()
    st.success('Done!')
    st.header(f'Amount of colonies detected: {check}')
    st.image(Image.fromarray(image_np))
