setwd("/_SONG/song-r/xc")

library(av)
library(tuneR)
library(seewave)
library(ggplot2)


# select the song cut from each recording (recordings can be directly downloaded using XC id numbers from Xeno-canto.org)
av_audio_convert('XC19010 - Blue-banded Pitta - Erythropitta arquata.mp3', 'XC19010 song.mp3', start_time = 4.6, total_time = 4.3)
av_audio_convert('XC451806 - Ash-throated Casiornis - Casiornis fuscus.mp3', 'XC451806 song.mp3', start_time = 15, total_time = 4.8)
av_audio_convert('XC100475 - Amazonian Antpitta - Hylopezus berlepschi.mp3', 'XC100475 song.mp3', start_time = 21.7, total_time = 4.8)
av_audio_convert('XC335063 - Yellow-breasted Antpitta - Grallaria flavotincta.mp3', 'XC335063 song.mp3', start_time = 9.58, total_time = 1.97)
av_audio_convert('XC275051 - Esmeraldas Antbird - Sipia nigricauda.mp3', 'XC275051 song.mp3', start_time = 5, total_time = 4.3)
av_audio_convert('XC229321 - Varzea Schiffornis - Schiffornis major.mp3', 'XC229321 song.mp3', start_time = 0.8, total_time = 6)

# read the song cuts
s1=readMP3('XC19010 song.mp3')
s2=readMP3('XC451806 song.mp3')
s3=readMP3('XC100475 song.mp3')
s4=readMP3('XC335063 song.mp3')
s5=readMP3('XC275051 song.mp3')
s6=readMP3('XC229321 song.mp3')

# set the colour palette
pal = spectro.colors(30)[c(10:1, 15:30)]

# draw the spectrograms (to show the change in frequency)
cus.spectrogram <- function(sound_file, A_limit, freq_limits, freq_breaks, xlab = '\n Time (s)', ylab = 'Frequency \n (kHz) \n '){
    ggspectro(sound_file, f=sound_file@samp.rate, wl=1024, wn="hanning", ovlp=87.5) + #
        geom_raster(aes(fill = amplitude), interpolate = T) +
        scale_fill_gradientn(colours = pal, limits=c(A_limit, 0), na.value = 'black') +
        theme_classic() + 
      scale_x_continuous(expand=c(0, 0)) + 
      scale_y_continuous(expand=c(0, 0), limits = freq_limits, breaks = freq_breaks) +
      theme(legend.position = 'none', 
            axis.line = element_blank(),
            axis.text = element_text(face="bold", size=12),
            axis.title = element_text(size = 14),
            plot.margin = margin(t=0, r=5, b=0, l=17)) +
      labs(x = xlab, y = ylab)
}

p1 = cus.spectrogram(s1, -36, c(0, 4.8), c(1:4), '')
p2 = cus.spectrogram(s2, -45, c(0, 7.9), c(2,4,6))
p3 = cus.spectrogram(s3, -42, c(0.01, 2.1), c(1,2), '', '')
p4 = cus.spectrogram(s4, -35, c(0, 5.8), c(2,4), ylab = '')
p5 = cus.spectrogram(s5, -43, c(1.5, 8.2), c(3,5,7), '', '')
p6 = cus.spectrogram(s6, -38, c(0, 5.8), c(2,4), ylab = '')



# draw the waveforms (to show the change in amplitude)
waveform = function(sound_file, ylab = '\nAmplitude \n(kU)', A){
    sample = seq(1:length(sound_file@left))
    time <-  sample / sound_file@samp.rate
    amplitude <- as.vector(cbind(sound_file@left))/1000
    df <- data.frame(sample, time, amplitude)
    
    ggplot(df)+
      geom_line(mapping = aes(x=time, y=amplitude), color="navy")+ 
      theme_classic() +
      scale_x_continuous(expand=c(0, 0)) + 
      scale_y_continuous(breaks = c(-A, 0 , A)) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x =element_blank(),
            axis.text.y = element_text(face="bold", size=12),
            axis.title = element_text(size = 14),
            axis.line = element_line(colour = 'grey'),
            plot.margin = margin(t=10, r=5, b=0, l=5)) +
      labs(x = '', y = ylab)
}

w1 = waveform(s1, A=20)
w2 = waveform(s2, A=20)
w3 = waveform(s3, '', A=10)
w4 = waveform(s4, '', A=10)
w5 = waveform(s5, '', A=10)
w6 = waveform(s6, '', A=10)


# add waveforms on top of the spectrograms
library(cowplot)
pdf('Fig.2.pdf', width = 14, height = 7)
plot_grid(w1, w3, w5, 
          p1, p3, p5, 
          w2, w4, w6, 
          p2, p4, p6, 
          ncol = 3, scale = 0.98,
          rel_heights = c(1, 2.4, 1, 2.6), rel_widths = c(1.1, 1, 1.05),
          labels = c('a', ' c', ' e', 
                     '', '', '',
                     'b', ' d', ' f'), label_size = 18)
dev.off()



