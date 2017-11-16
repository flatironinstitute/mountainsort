Waveform drift
============================

Spike sorting
-------------

Let $X$ be a $M\times N$ multi-channel timeseries array representing the signal acquired during a recording session. Here $M$ is the number of channels and $N$ is the number of timepoints in the recording. A spike sorting algorithm applied to this data will yield a sequence of firing events:


..  math::

  (t_1,k_1),\dots,(t_L,k_L)

where $t_j$ is the timepoint of the $j^{th}$ event and $k_j$ is the corresponding neuron label

..  math::

 k_j \in {1,\dots,K}.

Ideally, assuming that each firing event produces a signal spike whose shape depends only on neuron label, we would have

..  math::

  X(m,t)=\sum_{j=1}^{L} W_{k_j}(m,t-t_j) + \eta(m,t)


where $W_k$ is the waveform corresponding to the $k^{th}$ neuron, and $\eta$ is noise. We can estimate the average spike waveforms, or templates, from the output of spike sorting as follows:

..  math::

  W_k(m,t)=\sum_{j,k_j=k}{X(m,t_j+t)}

Drift
-----

In practice the spike waveform for a given neuron will change slightly over time. This is known as drift. Our signal model becomes more complicated:

..  math::

  X(m,t)=\sum_{j=1}^{L} W_{k_j,t_j}(m,t-t_j) + \eta(m,t)

where the waveform $W_{k,t}$ for the $k^{th}$ neuron now has a time dependency.

Images
------

.. image:: drift1.png

.. image:: drift2.png

.. image:: drift3.jpg

.. image:: drift4.jpg