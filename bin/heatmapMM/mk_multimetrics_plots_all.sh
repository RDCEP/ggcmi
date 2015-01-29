#!/bin/bash

gadm=(1 3 4 7 10 11 12 14 15 16 17 18 19 20 21 22 23 24 26 27 29 30 32 34 35 36 37 38 39 40 41 44 45 46 47 48 52 53 54 55 56 57 58 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 83 84 85 86 87 89 91 92 93 94 96 97 98 99 101 103 104 105 106 107 108 109 110 111 112 113 114 116 117 118 119 121 122 123 124 125 126 127 128 130 131 134 135 136 138 139 140 141 142 143 144 145 146 147 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 173 174 175 176 177 179 180 181 182 183 184 185 186 187 188 189 191 192 197 198 199 200 201 202 203 205 206 207 208 209 211 212 213 215 216 217 218 219 220 221 222 223 224 225 226 227 229 230 231 232 233 234 236 237 238 239 240 242 243 244 246 247 250 251 252 253)

for g in ${gadm[@]}; do
   for v in tscorr varratio rmse; do
      for m in none variance-scale mean-scale; do
         for c in mai ric soy whe mil sor; do
            echo $g, $v, $m, $c
            ./multimetrics_plot.py -d /project/ggcmi/AgMIP.output/processed/multimetrics/gadm0/faostat/fixed -v $v -x model -y climate,AgCFSR,AgMERRA,CFSR,ERAI,GRASP,Princeton,WATCH,WFDEI.CRU,WFDEI.GPCC -c gadm0,$g -o scen -c mp,true -c cr,$m -c time_range,full -c crop,$c -o dt --outdir /project/ggcmi/AgMIP.output/processed/plots/heatmapMM -f jpg
         done
      done
   done
done
