#!/usr/bin/env python3
import argparse
import glob
import json
import logging
import os
import re
import shlex
import subprocess
import sys
from collections import defaultdict

import plotly.graph_objects as go
from plotly.offline import plot
from plotly.subplots import make_subplots

import translate

logging.basicConfig(
    format='[%(asctime)s' '] %(levelname)s: %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M')

bin_dir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def setup_argparse():
    parser = argparse.ArgumentParser(description='''Script to do recombinant analysis''')
    parser.add_argument('-m', '--minAF',metavar='[FLOAT]',required=False, type=float, default=0.1, help="plot threshold. minimum average alleic frequency of strains' unqiue mutations [default:0.1]")
    parser.add_argument('--minRecAF',metavar='[FLOAT]',required=False, type=float, default=0.3, help="recobminant threshold. minimum average alleic frequency of strains' unqiue mutations [default:0.3]")
    parser.add_argument('-x', '--mixRatio',metavar='[FLOAT]',required=False, type=float, default=0.5, help="ratio of alleic frequency sites between 0.2 ((minMixAF)) and 0.8 (maxMixAF) to determine the mixed population")
    parser.add_argument('--minMixed_n',metavar='[INT]',required=False, type=int, default=3, help="threadhold of mixed mutations count.")
    parser.add_argument('--minMixAF',metavar='[FLOAT]',required=False, type=float, default=0.2, help="minimum alleic frequency for checking mixed mutation [default:0.2]")
    parser.add_argument('--maxMixAF',metavar='[FLOAT]',required=False, type=float, default=0.8, help="maximum alleic frequency for checking mixed mutation [default:0.8]")
    parser.add_argument('--lineageMutation',metavar='[FILE]',required=False, type=str, help="lineage mutation json file")
    parser.add_argument('--variantMutation',metavar='[FILE]',required=False, type=str, help="variant mutation json file")
    parser.add_argument('-eo', '--ec19_projdir',metavar='[PATH]',required=True, type=str,  help="ec-19 project directory")
    parser.add_argument('--igv', metavar='[PATH]',required=False, type=str,  help="igv.html relative path")
    parser.add_argument('--html' ,metavar='[FILE]',required=False, type=str, help='output plot html')
    parser.add_argument('--verbose', action='store_true',help='Show more infomration in log')
    argvs = parser.parse_args()

    if not argvs.lineageMutation:
        argvs.lineageMutation = os.path.join(bin_dir, 'lineage_mutation.json')
    if not argvs.variantMutation:
        argvs.variantMutation = os.path.join(bin_dir, 'variant_mutation.json')
    if not argvs.html:
        argvs.html = os.path.join(argvs.ec19_projdir, 'ReadsBasedAnalysis','readsMappingToRef','recombinant_analysis_result.html')
    if not argvs.igv:
        argvs.igv = os.path.join('..', '..','IGV','ref_tracks','igv.html')
    return argvs

def parse_variants(vcf,comp,nt_to_variant,argvs):
    mutations_list = list(nt_to_variant.keys())
    mutation_af=dict()
    mutation_dp=dict()
    pos_list = [u.split(':')[1] for u in mutations_list]
    comp_pos=dict()
    vcf_comp=dict()
    vcf_sep_comma=defaultdict(dict)
    parents_v=defaultdict(dict)
    mix_count_in_known_variants=dict()
    with open(comp,'r') as c:
        previous_content=[]
        for line in c:
            if line.startswith("#"):
                continue
            content = line.strip().split('\t')
            comp_pos[content[1]] = content

    vcf_content = []
    mix_count=0
    with open(vcf,'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            vcf_content.append(line)
            content = line.strip().split('\t')
            content2 = content[-1].split(':')
            AFreq = float(content2[-1])
            comp_content = comp_pos[content[1]]
            ref_bases = content[3]
            alt_bases = content[4]
            mut_list=[]
            if ',' in alt_bases or (AFreq > argvs.minMixAF and AFreq < argvs.maxMixAF):
                mix_count += 1
            if ',' in alt_bases:
                alt_list = alt_bases.split(',')
                if len(ref_bases) == 1:
                    for alt in alt_list:
                        if len(alt) == 1:
                            mut = ref_bases + ":" + str(int(content[1])) + ":" + str(alt)
                        if len(alt) > 1:
                            comp_content = comp_pos[str(int(content[1]) + 1)]
                            mut = "ins" + ":" + str(int(content[1]) + 1) + ":" + str(alt[1:])
                        mut_list.append(mut)
                else: # len(ref_bases) > 1:
                    for alt in alt_list:
                        if len(alt) < len(ref_bases):
                            comp_content = comp_pos[str(int(content[1]) + 1)]
                            mut = 'del' + ":" + str(int(content[1]) + 1) + ":" + str(len(ref_bases)-len(alt))
                        if len(alt) > len(ref_bases):
                            comp_content = comp_pos[str(int(content[1]) + 1)]
                            mut = "ins" + ":" + str(int(content[1]) + 1) + ":" + str(alt.replace(ref_bases,''))
                        if len(alt) == len(ref_bases):
                            if alt[0] == ref_bases[0]:
                                alt=alt[1:]
                                ref_bases=ref_bases[1:]
                            mut = ref_bases + ":" + str(int(content[1]) + 1) + ":" + str(alt)
                        mut_list.append(mut)
                vcf_comp[ref_bases + ":" + content[1] + ":" + alt_bases] = ':'.join([ comp_content[x].split(' ')[0] for x in range(4,9) ])
            else:
                if len(ref_bases) - len(alt_bases) != 0:
                    comp_content = comp_pos[str(int(content[1]) + 1)]
                    ref_bases = 'del' if len(content[3]) - len(content[4]) > 0 else 'ins'
                    alt_bases = abs(len(content[3]) - len(content[4])) if ref_bases == 'del' else content[4].replace(content[3],'')
                    mut = ref_bases + ":" + str(int(content[1]) + 1) + ":" + str(alt_bases)
                    vcf_comp[mut] = ':'.join([ comp_content[x].split(' ')[0] for x in range(4,9) ])    
                else:
                    mut = ref_bases + ":" + str(int(content[1])) + ":" + str(alt_bases)
                    vcf_comp[mut] = ':'.join([ comp_content[x].split(' ')[0] for x in range(4,9) ])
                
                mut_list.append(mut)
            for nt_v in mut_list:
                vcf_sep_comma[nt_v]['AF']=AFreq
                vcf_sep_comma[nt_v]['DP']=comp_content[3]
                if nt_v in nt_to_variant:
                    vcf_sep_comma[nt_v]['variants']=nt_to_variant[nt_v]
                    vcf_sep_comma[nt_v]['DP']=comp_content[3]
                    if ',' in str(alt_bases) or (AFreq > argvs.minMixAF and AFreq < argvs.maxMixAF):
                        mix_count_in_known_variants[nt_v] = 1
                    ## scan first time find unique mutations variants
                    if len(nt_to_variant[nt_v]) == 1:
                        if ''.join(nt_to_variant[nt_v]) in parents_v:
                            parents_v[''.join(nt_to_variant[nt_v])]['uniq'] += 1 
                        else:
                            parents_v[''.join(nt_to_variant[nt_v])]['uniq'] = 1
                            parents_v[''.join(nt_to_variant[nt_v])]['all'] = 0
                            parents_v[''.join(nt_to_variant[nt_v])]['filtered'] = 0
                else:
                    vcf_sep_comma[nt_v]['variants']=[]
    # use the first scan, uniq list to count all mutataions with the variants and variants between AF range
    for v in parents_v:
        for nt_v in vcf_sep_comma:
            if v in vcf_sep_comma[nt_v]['variants']:
                parents_v[v]['all'] += 1
                if vcf_sep_comma[nt_v]['AF'] > argvs.minMixAF and vcf_sep_comma[nt_v]['AF'] < argvs.maxMixAF:
                    parents_v[v]['filtered'] += 1
    parents_v = dict(sorted(parents_v.items(), key=lambda item: item[1]['all'], reverse=True))

    for u in mutations_list:
        exist = 0
        ref_nt,pos,alt_nt = u.split(':')
   
        for v in vcf_content:
            content = v.strip().split('\t')
            content2 = content[-1].split(':')
            AFreq = content2[-1]
            #total_dp = int(content2[1]) + int(content2[2])
            if ref_nt == 'del' and str(int(pos) - 1) == content[1] and str(len(content[3]) - len(content[4])) == alt_nt:
                #logging.info( "\t".join(content) +  str(len(content[3]) - len(content[4])) + "\t" + alt_nt + "\n" )               
                mutation_af[u] = AFreq
                exist = 1 
            elif ref_nt == 'ins' and str(int(pos) - 1) == content[1] and alt_nt in content[4]:
                mutation_af[u] = AFreq
                exist = 1
            elif pos == content[1] and ref_nt == content[3] and alt_nt == content[4]:
                mutation_af[u] = AFreq
                exist = 1
        comp_content = comp_pos[pos]
        all_count = ':'.join([ comp_content[x].split(' ')[0] for x in range(4,9) ])
        mutation_dp[u]=str(all_count)
        if not exist:
            (count, percentage) = (str(0), str(0))
            total_dp = comp_content[3]
            if ref_nt == 'del' or ref_nt == 'ins':
                count, percentage=comp_content[8].split(' ')
            else:
                if alt_nt == 'A':
                    count, percentage=comp_content[4].split(' ')
                if alt_nt == 'C':
                    count, percentage=comp_content[5].split(' ')
                if alt_nt == 'G':
                    count, percentage=comp_content[6].split(' ')
                if alt_nt == 'T':
                    count, percentage=comp_content[7].split(' ')
                if not alt_nt:
                    ref_base = comp_content[2]
                    if ref_base == 'A':
                        count, percentage=comp_content[4].split(' ')
                    if ref_base == 'C':
                        count, percentage=comp_content[5].split(' ')
                    if ref_base == 'G':
                        count, percentage=comp_content[6].split(' ')
                    if ref_base == 'T':
                        count, percentage=comp_content[7].split(' ')
                    percentage = "{:.4f}".format( 1 - (int(count)/int(total_dp))) if int(total_dp) > 0 else '0.000' 

            mutation_af[u]=percentage.replace('(','').replace(')','')
        if u in vcf_sep_comma:
            vcf_sep_comma[u]['AF'] = mutation_af[u]

    mutations_count = len(vcf_comp)
    mix_ratio = mix_count/len(vcf_comp) * 100 if len(vcf_comp) > 0 else 0
    logging.info(f"{mix_count}/{mutations_count} mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF}.")
    logging.info(f"{len(mix_count_in_known_variants)}/{mix_count} mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF} and known variants.")
    logging.info(f"All probable parents, mutation count: {parents_v}")
    if len(vcf_comp) > 0 and mix_count/len(vcf_comp) > argvs.mixRatio:
        logging.info(f"Probable Mixed Infection. {mix_ratio:.2f}% ({mix_count}/{mutations_count}) mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF}.")

    return mutation_af, mutation_dp, vcf_comp, mix_count, vcf_sep_comma, parents_v

def check_mutations(mutation_af,mutation_dp, delta_uniq_nt,omicron_uniq_nt,argvs):
    ## need to set criteria to plot.  default to False.
    check=False
    probably_delta = 0
    probably_omicron = 0
    total_delta_af = 0
    total_omicron_af = 0
    for u in delta_uniq_nt:
        ac,cc,gc,tc,nc = mutation_dp[u].split(':')
        total_delta_af = total_delta_af + float(mutation_af[u])
        if float(mutation_af[u]) >= argvs.maxMixAF:
            probably_delta += 1
        if float(mutation_af[u]) <= argvs.minMixAF and (int(ac) + int(cc) + int(gc) + int(tc) + int(nc)) > 0:
            probably_omicron += 1
    for u in omicron_uniq_nt:
        ac,cc,gc,tc,nc = mutation_dp[u].split(':')
        total_omicron_af = total_omicron_af + float(mutation_af[u])
        if float(mutation_af[u]) >= argvs.maxMixAF:
            probably_omicron += 1
        if float(mutation_af[u]) <= argvs.minMixAF and (int(ac) + int(cc) + int(gc) + int(tc) + int(nc)) > 0:
            probably_delta += 1
    avg_delta_AF = total_delta_af/len(delta_uniq_nt)
    avg_omicron_AF = total_omicron_af/len(omicron_uniq_nt)
    #if probably_delta > float(minAF) * len(delta_uniq_nt) or probably_omicron > float(minAF) * len(omicron_uniq_nt):
    if avg_delta_AF > argvs.minAF or avg_omicron_AF > argvs.minAF:
        check=True
    logging.info(f'Average Delta Unique Mutations AF: {avg_delta_AF}')
    logging.info(f'Average Omicron Unique Mutations AF: {avg_omicron_AF}')
    logging.info(f'Probably Delta Unique Mutations Count (>={argvs.maxMixAF}): {probably_delta}')
    logging.info(f'Probably Omicron Unique Mutations Count (>={argvs.maxMixAF}): {probably_omicron}')
    if avg_delta_AF > argvs.minRecAF and avg_omicron_AF > argvs.minRecAF and int(len(delta_uniq_nt) * argvs.minRecAF) < probably_delta  and int(len(omicron_uniq_nt) * argvs.minRecAF) < probably_omicron :
        logging.info(f'Probable Delta and Omicron recombinant. [minimum average delta and omicron unqiue variants AF: {argvs.minRecAF}]')
    return check

def load_var_mutation(file):
    with open( file, 'r') as f:
        v_data = json.load(f)
    delta_uniq_nt=list(set(v_data['Delta'].keys()) - set(v_data['Omicron'].keys()))
    omicron_uniq_nt=list(set(v_data['Omicron'].keys()) - set(v_data['Delta'].keys()))
    
    nt_to_aa=dict()
    nt_to_variant=dict()
    for k in v_data:
        nt_to_aa.update(v_data[k])
        for nt in v_data[k]:
            if nt not in nt_to_variant:
                nt_to_variant[nt] = [k]
            else:
                nt_to_variant[nt].append(k)

    return(delta_uniq_nt,omicron_uniq_nt,nt_to_variant,nt_to_aa)

def load_lineage_mutation(file):
    with open( file, 'r') as f:
        l_data = json.load(f)
    nt_to_lineage=dict()
    nt_to_aa=dict()
    for k in l_data:
        nt_to_aa.update(l_data[k])
        for nt in l_data[k]:
            if nt not in nt_to_lineage:
                nt_to_lineage[nt] = [k]
            else:
                nt_to_lineage[nt].append(k)
    return nt_to_lineage

def deltacron_af_plot(nt_to_aa, delta_uniq_nt,omicron_uniq_nt,mutation_af,barplot_lists, lineage , sample, url , output):
    new_d = { i: int(i.split(':')[1]) for i in delta_uniq_nt + omicron_uniq_nt }
    deltacron_variants_nt = list(dict(sorted(new_d.items(), key=lambda item: item[1])).keys())
    deltacron_variants_aa = [ nt_to_aa[x] for x in deltacron_variants_nt]
    fig = make_subplots(rows=2,cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[200,500])
   
    fig.add_trace(go.Scatter(
        x= deltacron_variants_nt,
        y=[ float(mutation_af[i]) if i in delta_uniq_nt else None for i in deltacron_variants_nt ],
        mode="markers",
        marker=dict(color='red',size=8),
        name="Delta unique variants",
        hovertemplate = 'AF: %{y:.2f}<br>'+'<b>%{text}</b>',
        text=[ f'{i}' for i in deltacron_variants_aa ],
        customdata=[ url + '?locus=NC_045512_2:' + str(int(i.split(':')[1]) - 100) + '-' +  str(int(i.split(':')[1]) + 100) if i in delta_uniq_nt else None for i in deltacron_variants_nt ],
        showlegend=True,
        legendgroup=1
    ),row=2,col=1)
    fig.add_trace(go.Scatter(
        x=deltacron_variants_nt,
        y=[ float(mutation_af[i]) if i in omicron_uniq_nt else None for i in deltacron_variants_nt ],
        mode="markers",
        marker=dict(color='blue',size=8),
        name="Omicron unique variants",
        hovertemplate = 'AF: %{y:.2f}<br>'+'<b>%{text}</b>',
        text=[ f'{i}' for i in deltacron_variants_aa ],
        customdata=[ url + '?locus=NC_045512_2:' + str(int(i.split(':')[1]) - 100) + '-' +  str(int(i.split(':')[1]) + 100) if i in omicron_uniq_nt else None for i in deltacron_variants_nt ],
        showlegend=True,
        legendgroup=1
    ),row=2,col=1)
    fig.add_trace(go.Bar(
        x=deltacron_variants_nt,
        y=barplot_lists[0],
        marker_color='red',
        name="Delta Depth Coverage",
        showlegend=False,
        offsetgroup=0,
    ),row=1,col=1)
    fig.add_trace(go.Bar(
        x=deltacron_variants_nt,
        y=barplot_lists[1],
        marker_color='blue',
        name="Omicron Depth Coverage",
        showlegend=False,
        offsetgroup=0,
        base=barplot_lists[0]
    ),row=1,col=1)
    fig.add_trace(go.Bar(
        x=deltacron_variants_nt,
        y=barplot_lists[2],
        marker_color='grey',
        name="Other Depth Coverage",
        showlegend=False,
        offsetgroup=0,
        base=[val1+val2 for val1, val2 in zip(barplot_lists[0],barplot_lists[1])],
    ),row=1,col=1)
    fig.add_vline(x=10.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=14.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=42.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=43.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=44.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=48.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=50.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=51.5, line_width=1, line_dash="dash", line_color="grey")
            
    ## update x tick val to amino acid change
    #fig.update_layout(xaxis=dict(tickvals=vendor_names,ticktext=[ '%s (n=%d)'% (vendor_names[m], vendor_count[m]) for m in range(len(vendor_count))] ),font=dict(size=16))
    fig.update_yaxes(range=[-0.1, 1.1],title='Alternative Frequency',row=2,col=1)
    fig.update_yaxes(title='Depth Coverage',row=1,col=1)
    #fig.update_layout(barmode='stack',row=1,col=1)
    #fig.update_xaxes(tickvals = list(range(len(variants))), ticktext=variants_aa, tickfont=dict(size=8),tickangle=45)
    fig.update_layout(
        title=sample + ' (' + lineage + ')',
        #height=1800, width=1200,
        #width=1400,
        legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.0,
                xanchor="right",
                x=0.95
            ),
    )
    # Get HTML representation of plotly.js and this figure
    plot_div = plot(fig, output_type='div', include_plotlyjs=True)

    # Get id of html div element that looks like
    # <div id="301d22ab-bfba-4621-8f5d-dc4fd855bb33" ... >
    res = re.search('<div id="([^"]*)"', plot_div)
    div_id = res.groups()[0]

    # Build JavaScript callback for handling clicks
    # and opening the URL in the trace's customdata 
    js_callback = """
    <script>
    var plot_element = document.getElementById("{div_id}");
    plot_element.on('plotly_click', function(data){{
        var point = data.points[0];
        if (point) {{
            window.open(point.customdata);
        }}
    }})
    </script>
    """.format(div_id=div_id)

    # Build HTML string
    html_str = """
    <html>
    <body>
    {plot_div}
    {js_callback}
    </body>
    </html>
    """.format(plot_div=plot_div, js_callback=js_callback)
    with open(output, 'w') as f:
        f.write(html_str)
    fig.write_image(output+'.png')
    return output

def bar_plot_list3(delta_uniq_nt,omicron_uniq_nt,mutation_dp):
    count_dict=dict()
    delta_count_list=list()
    omicron_count_list=list()
    other_count_list=list()
    new_d = { i: int(i.split(':')[1]) for i in delta_uniq_nt + omicron_uniq_nt }
    deltacron_variants_nt = list(dict(sorted(new_d.items(), key=lambda item: item[1])).keys())
    for i in deltacron_variants_nt:
        ref_nt,pos,alt_nt = i.split(':')
        ac,cc,gc,tc,nc = mutation_dp[i].split(':')
        count_dict['A']=int(ac)
        count_dict['C']=int(cc)
        count_dict['G']=int(gc)
        count_dict['T']=int(tc)
        count_dict['N']=int(nc)
        total_nt = int(ac) + int(cc) + int(gc) + int(tc) + int(nc)
        if i in delta_uniq_nt:
            if pos == '23604':
                delta_count = count_dict[alt_nt]
                omicron_count = count_dict['A']
                other_count = total_nt - delta_count - omicron_count
            elif pos == '28881':
                delta_count = count_dict[alt_nt]
                omicron_count = count_dict['A']
                other_count = total_nt - delta_count - omicron_count
            elif alt_nt in ['A','T','G','C']:
                delta_count = count_dict[alt_nt]
                omicron_count = count_dict[ref_nt]
                other_count = total_nt - delta_count - omicron_count
            else:
                delta_count = count_dict['N']
                omicron_count = total_nt - delta_count
                other_count = 0
        elif i in omicron_uniq_nt:
            if pos == '23604':
                omicron_count = count_dict[alt_nt]
                delta_count = count_dict['G']
                other_count = total_nt - delta_count - omicron_count
            elif pos == '28881':
                omicron_count = count_dict[alt_nt]
                delta_count = count_dict['T']
                other_count = total_nt - delta_count - omicron_count
            elif alt_nt in ['A','T','G','C']:
                delta_count = count_dict[ref_nt]
                omicron_count = count_dict[alt_nt]
                other_count = total_nt - delta_count - omicron_count
            else:
                omicron_count = count_dict['N']
                delta_count = total_nt - omicron_count
                other_count = 0
        else:
            omicron_count = None
            delta_count = None
            other_count = None
        delta_count_list.append(delta_count)
        omicron_count_list.append(omicron_count)
        other_count_list.append(other_count)
    return [delta_count_list,omicron_count_list,other_count_list]

def deltacron_af_plot_by_sample_id(nt_to_variant, nt_to_aa, delta_uniq_nt,omicron_uniq_nt, sample, mutation_af, mutation_dp, lineage, argvs):
    igv_url= argvs.igv
    out_html = os.path.splitext(argvs.html)[0] + '.deltacron.html'
    logging.info(f"Generating deltacron unique mutations AF plot and save to {out_html}")
    barplot_lists = bar_plot_list3(delta_uniq_nt,omicron_uniq_nt, mutation_dp)
    deltacron_af_plot(nt_to_aa,delta_uniq_nt,omicron_uniq_nt,mutation_af, barplot_lists, lineage, sample, igv_url,out_html)

def comp_stack_bar_plot(vcf_comp, mix_count, nt_to_lineage, projname, lineage, argvs):
    barplot_html_output = os.path.join(os.path.dirname(argvs.html),'variants_bar_plot.html')
    logging.info(f"Generating mutations bar plot and save to {barplot_html_output}")
    count_dict=dict()
    #vcf_list = list(vcf_comp.keys())
    igvurls =[ argvs.igv + '?locus=NC_045512_2:' + str(int(i.split(':')[1]) - 100) + '-' +  str(int(i.split(':')[1]) + 100) for i in vcf_comp.keys() ]
    vcf_list = [ i.split(":")[0] + ":<b>" + i.split(":")[1] + "</b>:" + i.split(":")[2] for i in vcf_comp.keys()]
    aclist, cclist, gclist, tclist, nclist = list(), list(),list(),list(),list()
    aplist, cplist, gplist, tplist, nplist = list(), list(),list(),list(),list()
    mutations_count = len(vcf_comp)
    mix_ratio = mix_count/mutations_count * 100 if mutations_count > 0 else 0
    variants_list = []
    variants_list_val = []
    variants_list_color = []
    variants_list_labels = []
    for i in vcf_comp:
        ac,cc,gc,tc,nc = vcf_comp[i].split(':')
        total = int(ac) + int(cc) + int(gc) + int(tc) + int(nc)
        aclist.append(int(ac))
        aplist.append(int(ac)/total*100)
        cclist.append(int(cc))
        cplist.append(int(cc)/total*100)
        gclist.append(int(gc))
        gplist.append(int(gc)/total*100)
        tclist.append(int(tc))
        tplist.append(int(tc)/total*100)
        nclist.append(int(nc))
        nplist.append(int(nc)/total*100)
        mut_list=[i,None]
        if ',' in i:
            mut_list = i.split(',')
            if len(mut_list[0].split(":")[0]) == 1:
                mut_list[1] = ':'.join(mut_list[0].split(":")[0:2] + [mut_list[1]])
            if len(mut_list[0].split(":")[0]) == 2:
                if len(mut_list[0].split(":")[2]) == 1:
                    mut_list[0] =  "del:" + str(int(mut_list[0].split(":")[1]) + 1) + ":1"
                if len(mut_list[0].split(":")[2]) == 2:
                    mut_list[0] =  mut_list[0].split(":")[0][1] + ":" + str(int(mut_list[0].split(":")[1]) + 1) + ":" + mut_list[0].split(":")[2][1]
                if len(mut_list[1]) == 1:
                    mut_list[1] =  "del:" + str(int(mut_list[0].split(":")[1]) + 1) + ":1" 
                if len(mut_list[1]) == 2:
                    mut_list[1] =  mut_list[0].split(":")[0][1] + ":" + str(int(mut_list[0].split(":")[1]) + 1) + ":" + mut_list[1][1]

            #print(mut_list)
        if mut_list[0] in nt_to_lineage or mut_list[1] in nt_to_lineage:
            variants_list.append(1)
            variants_list_val_text_dict=dict()
            variants_list_val_text=""
            variants_list_label_text=""
            red_n=0
            green_n=0
            blue_n=0
            var_list=[]
            if mut_list[0] in nt_to_lineage:
                var_list.extend(nt_to_lineage[mut_list[0]])
            if mut_list[1] in nt_to_lineage:
                var_list.extend(nt_to_lineage[mut_list[1]])
            for vi, var in enumerate(var_list):
                if var.startswith('B.1.1.529') or var.startswith('BA'):
                    #Omicron
                    variants_list_val_text_dict[var]=1
                    blue_n=255
                elif var.startswith('AY') or var.startswith('B.1.617.2'):
                    # Delta
                    var1 = re.sub(r"(AY\.\d+)\.\d+\.?\d?",r"\1",var)
                    variants_list_val_text_dict[var1]=1
                    red_n=255
                elif var == 'B.1.1.7' or var.startswith('Q.') or var.startswith('B.1.351') or var=='P.1' or var.startswith('B.1.427') or var.startswith('B.1.429') or var.startswith('B.1.525') or var.startswith('B.1.526') or var.startswith('B.1.671') or var.startswith('B.1.621') or var == 'P.2':
                    # alpha(B1.1.7 and Q.*), beta (B.1.351.*), gamma (P.1), epsilon, Eta, lota, kappa, mu, zeta
                    var1 = re.sub(r"(P\.1)\.\d+\.?\d?",r"\1",var)
                    var2 = re.sub(r"(B\.1\.351)\.?\d?",r"\1",var1)
                    var3 = re.sub(r"(B\.1\.621)\.?\d?",r"\1",var2)
                    var4 = re.sub(r"(Q.\d+)",r"Q.x",var3)
                    variants_list_val_text_dict[var4]=1
                    green_n=180
                else:
                    # other sublineage
                    variants_list_val_text_dict['Others']=1
                    green_n=180

            for vi, var in  enumerate(sorted(list(variants_list_val_text_dict.keys()))):
                if vi % 3 == 2:
                    variants_list_val_text += var + "<br>"
                else:
                    variants_list_val_text += var + " ,"        
            variants_list_val.append(variants_list_val_text)                
            variants_list_color.append(f"rgba({red_n},{green_n},{blue_n},1)")
            if blue_n == 255:
                variants_list_label_text += 'O'
            if red_n == 255:
                variants_list_label_text += 'D'
            if green_n == 180:
                variants_list_label_text += 'N'
            variants_list_labels.append(variants_list_label_text)
        else:
            variants_list.append(None)
            variants_list_val.append(None)
            variants_list_labels.append(None)
            variants_list_color.append('rgba(0,0,0,0)')

    fig = make_subplots(rows=3,cols=1, shared_xaxes=True, vertical_spacing=0.01, row_heights=[200,500,15])
    fig.add_trace(go.Bar(name='A', x=vcf_list, y=aclist,marker_color='green',opacity=.7),row=1,col=1)
    fig.add_trace(go.Bar(name='C', x=vcf_list, y=cclist,marker_color='blue',opacity=.7),row=1,col=1)
    fig.add_trace(go.Bar(name='G', x=vcf_list, y=gclist,marker_color='orange',opacity=.7),row=1,col=1)
    fig.add_trace(go.Bar(name='T', x=vcf_list, y=tclist,marker_color='red',opacity=.7),row=1,col=1)
    fig.add_trace(go.Bar(name='INDELs/N', x=vcf_list, y=nclist,marker_color='grey',opacity=.7),row=1,col=1)

    fig.add_trace(go.Bar(name='A', x=vcf_list, y=aplist,marker_color='green',showlegend=False,opacity=.7, text=[ f"{i:.2f}" for i in aplist]),row=2,col=1)
    fig.add_trace(go.Bar(name='C', x=vcf_list, y=cplist,marker_color='blue', showlegend=False,opacity=.7, text=[ f"{i:.2f}" for i in cplist]),row=2,col=1)
    fig.add_trace(go.Bar(name='G', x=vcf_list, y=gplist,marker_color='orange', showlegend=False,opacity=.7, text=[ f"{i:.2f}" for i in gplist]),row=2,col=1)
    fig.add_trace(go.Bar(name='T', x=vcf_list, y=tplist,marker_color='red', showlegend=False,opacity=.7, text=[ f"{i:.2f}" for i in tplist]),row=2,col=1)
    fig.add_trace(go.Bar(name='INDELs/N', x=vcf_list, y=nplist,marker_color='grey', showlegend=False,opacity=.7, text=[ f"{i:.2f}" for i in nplist]),row=2,col=1)

    fig.add_trace(go.Bar(name='Variants', x=vcf_list, y=variants_list, marker_color=variants_list_color, text=variants_list_labels,textfont_size=6, textangle=-90,
                         hoverinfo='text', hoverlabel=dict(font_size=8), hovertext=variants_list_val,showlegend=False,opacity=1,
                         customdata=igvurls,
                        ),row=3,col=1)
  

    fig.update_layout(barmode='stack', title_text='Positions with mutations' + '<br><sup>' + projname + ' (' + lineage + ')</sup>')
    fig.update_xaxes(tickfont=dict(size=8),tickangle=-60)
    fig.update_yaxes(title="D.P for each nucleotide/indels",row=1,col=1)
    fig.update_yaxes(title="Percentage",row=2,col=1)
    fig.update_yaxes(visible=False,row=3,col=1)
    if mix_count > 0:
        fig.add_annotation(text=f"{mix_ratio:.2f}% ({mix_count}/{mutations_count}) mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF}.",
                  xref="x domain", yref="y domain", showarrow=False, font=dict(size=10),
                  x=0, y=1.12,row=1,col=1)
    #fig.write_html(barplot_html_output)
    fig.write_image(barplot_html_output+'.png')   
     # Get HTML representation of plotly.js and this figure
    plot_div = plot(fig, output_type='div', include_plotlyjs=True)

    # Get id of html div element that looks like
    # <div id="301d22ab-bfba-4621-8f5d-dc4fd855bb33" ... >
    res = re.search('<div id="([^"]*)"', plot_div)
    div_id = res.groups()[0]

    # Build JavaScript callback for handling clicks
    # and opening the URL in the trace's customdata 
    js_callback = """
    <script>
    var plot_element = document.getElementById("{div_id}");
    plot_element.on('plotly_click', function(data){{
        var point = data.points[0];
        if (point) {{
            window.open(point.customdata);
        }}
    }})
    </script>
    """.format(div_id=div_id)

    # Build HTML string
    html_str = """
    <html>
    <body>
    {plot_div}
    {js_callback}
    </body>
    </html>
    """.format(plot_div=plot_div, js_callback=js_callback)
    with open(barplot_html_output, 'w') as f:
        f.write(html_str)   

def mutations_af_plot(parents_v,vcf_sep_comma,nt_to_aa_class,projname,lineage, argvs):
    output = argvs.html
    logging.info(f"Generating mutations AF plots and save to {output}")
    all_mut_nt = [ i.split(":")[0] + ":<b>" + i.split(":")[1] + "</b>:" + i.split(":")[2] for i in list(vcf_sep_comma.keys())]
    
    parents = list(parents_v.keys())[0:2]
    fig = go.Figure()
    color1='blue'
    color2='red' 
    if parents[0] == 'Omicron' or parents[1] == 'Delta':
        color1='blue'
        color2='red'
    if parents[1] == 'Omicron' or parents[0] == 'Delta':
        color1='red'
        color2='blue'
    fig.add_trace(go.Scatter(
                x=all_mut_nt, 
                y=[ float(vcf_sep_comma[x]['AF'])  if parents[0] in vcf_sep_comma[x]['variants'] and parents[1] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma], 
                mode='markers',
                marker=dict(color=color1,size=8),
                name=parents[0],
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ vcf_sep_comma[x]['DP'] + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if parents[0] in vcf_sep_comma[x]['variants'] and parents[1] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma],
                customdata=[ argvs.igv + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if parents[0] in vcf_sep_comma[x]['variants'] and parents[1] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma ],
                showlegend=True,
                ))
    fig.add_trace(go.Scatter(
                x=all_mut_nt, 
                y=[ float(vcf_sep_comma[x]['AF'])  if parents[1] in vcf_sep_comma[x]['variants'] and parents[0] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma], 
                mode='markers',
                marker=dict(color=color2,size=8),
                name=parents[1],
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ vcf_sep_comma[x]['DP'] + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if parents[1] in vcf_sep_comma[x]['variants'] and parents[0] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma],
                customdata=[ argvs.igv + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if parents[1] in vcf_sep_comma[x]['variants'] and parents[0] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma ],
                showlegend=True,
                ))
    fig.add_trace(go.Scatter(
                x=all_mut_nt, 
                y=[ float(vcf_sep_comma[x]['AF'])  if parents[0] in vcf_sep_comma[x]['variants'] and parents[1] in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma], 
                mode='markers',
                marker=dict(color='purple',size=8),
                name=parents[0] + ', ' + parents[1],
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ vcf_sep_comma[x]['DP'] + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if parents[0] in vcf_sep_comma[x]['variants'] and parents[1] in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma],
                customdata=[ argvs.igv + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if parents[0] in vcf_sep_comma[x]['variants'] and parents[1] in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma ],
                showlegend=True,
                ))
    fig.add_trace(go.Scatter(
                x=all_mut_nt, 
                y=[ float(vcf_sep_comma[x]['AF'])  if vcf_sep_comma[x]['variants'] and parents[0] not in vcf_sep_comma[x]['variants'] and parents[1] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma], 
                mode='markers',
                marker=dict(color='green',size=8),
                name='Not ' + parents[0] + ', Not ' + parents[1],
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ vcf_sep_comma[x]['DP'] + '<br>' + 'AA: ' + nt_to_aa_class.convert_nt_prot(x) + '<br>Var: ' + ','.join(vcf_sep_comma[x]['variants']) if vcf_sep_comma[x]['variants'] and  parents[0] not in vcf_sep_comma[x]['variants'] and parents[1] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma],
                customdata=[ argvs.igv + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if vcf_sep_comma[x]['variants'] and parents[0] not in vcf_sep_comma[x]['variants'] and parents[1] not in vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma ],
                showlegend=True,
                ))
    fig.add_trace(go.Scatter(
                x=all_mut_nt, 
                y=[ float(vcf_sep_comma[x]['AF'])  if not vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma], 
                mode='markers',
                marker=dict(color='grey',size=8),
                name='Undefined Mutations',
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ vcf_sep_comma[x]['DP'] + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if not vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma],
                customdata=[ argvs.igv + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if not vcf_sep_comma[x]['variants'] else None for x in vcf_sep_comma ],
                showlegend=True,
                ))
    fig.update_xaxes(tickfont=dict(size=8),tickangle=-60)
    fig.update_yaxes(range=[-0.1, 1.1],title='Alternative Frequency')
    fig.update_layout(title="Positions with mutations" + '<br><sup>' + projname + ' (' + lineage + ')</sup>')
    # Get HTML representation of plotly.js and this figure
    plot_div = plot(fig, output_type='div', include_plotlyjs=True)

    # Get id of html div element that looks like
    # <div id="301d22ab-bfba-4621-8f5d-dc4fd855bb33" ... >
    res = re.search('<div id="([^"]*)"', plot_div)
    div_id = res.groups()[0]

    # Build JavaScript callback for handling clicks
    # and opening the URL in the trace's customdata 
    js_callback = """
    <script>
    var plot_element = document.getElementById("{div_id}");
    plot_element.on('plotly_click', function(data){{
        var point = data.points[0];
        if (point) {{
            window.open(point.customdata);
        }}
    }})
    </script>
    """.format(div_id=div_id)

    # Build HTML string
    html_str = """
    <html>
    <body>
    {plot_div}
    {js_callback}
    </body>
    </html>
    """.format(plot_div=plot_div, js_callback=js_callback)
    with open(output, 'w') as f:
        f.write(html_str)
    fig.write_image(output+'.png')
    return output

def parse_lineage(txt):
    with open(txt,'r') as f:
        content = f.readlines()
    result = content[1].split(',')
    return(result[1])

def parse_ec19_config(config):
    ec19_config=dict()
    with open(config,'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.__contains__('='):
                k,v =line.strip().split('=')
                ec19_config[k]=v
    return(ec19_config)

def process_cmd(cmd, msg=None, stdout_log=None,  shell_bool=False):
    if msg:
        logging.info(msg)

    logging.debug("CMD: %s" %(" ".join(cmd)) )
    proc = subprocess.Popen(shlex.split(" ".join(cmd)), shell=shell_bool, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    
    if stdout_log:
        f=open (stdout_log, "a") 
    while True:
        output = proc.stdout.readline()
        if output == '' and proc.poll() is not None:
            break
        if output:
            stdout_msg=output.decode()
            sys.stdout.write(stdout_msg)
            if stdout_log:
                f.write(stdout_msg)

    f.close()
    rc = proc.poll()
    if rc != 0: 
        logging.error("Failed %d %s" % (rc, " ".join(cmd)))
    
    return 0

def main():
    argvs = setup_argparse()
    log_level = logging.DEBUG if argvs.verbose else logging.INFO
    logging.basicConfig(
    format='[%(asctime)s' '] %(levelname)s: %(message)s', level=log_level, datefmt='%Y-%m-%d %H:%M')
    (delta_uniq_nt, omicron_uniq_nt, nt_to_variant, nt_to_aa) = load_var_mutation(argvs.variantMutation)
    nt_to_lineage = load_lineage_mutation(argvs.lineageMutation)
    nt_to_aa_class = translate.NT_AA_POSITION_INTERCHANGE(os.path.join(bin_dir, 'SARS-CoV-2.json'))

    ec19_consensus = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.fasta",recursive = False)
    ec19_log = glob.glob(argvs.ec19_projdir+f"/process.log",recursive = False)
    ec19_config_file = glob.glob(argvs.ec19_projdir+f"/config.txt",recursive = False)
    ec19_lineage_file = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus_lineage.txt",recursive = False)
    ec19_vcf = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.vcf",recursive = False)
    ec19_bam = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2.sort.bam",recursive = False)
    ec19_snp_result = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.SNPs_report.txt",recursive = False)
    ec19_indel_result = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.Indles_report.txt",recursive = False)
    ec19_compositionlog = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.compositionlog",recursive = False)
    ec19_aln_stats = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2.alnstats.txt",recursive = False)
    ec19_fq_count = glob.glob(argvs.ec19_projdir+f"/QcReads/fastqCount.txt",recursive = False)
    ec19_lineage_abund = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/LineageAbundance/predicitions.tsv",recursive = False)
    
    if not ec19_vcf:
        logging.error("Cannot find consensus.vcf file\n")
        sys.exit(1)
    if not ec19_compositionlog:
        logging.error("Cannot find consensus.compositionlog\n")
        sys.exit(1)

    ec19_lineage = parse_lineage(ec19_lineage_file[0]) if ec19_lineage_file else 'Unknown'
    ec19_config = parse_ec19_config(ec19_config_file[0])
    mutations_af, mutations_dp, vcf_comp, mix_count,vcf_sep_comma, parents_v = parse_variants(ec19_vcf[0],ec19_compositionlog[0],nt_to_variant,argvs)

    comp_stack_bar_plot(vcf_comp,mix_count,nt_to_lineage, ec19_config['projname'], ec19_lineage, argvs)

    #plot_bool = check_mutations(mutations_af, mutations_dp, delta_uniq_nt,omicron_uniq_nt,argvs)

    if mix_count > argvs.minMixed_n and len(parents_v.keys()) > 1:
        mutations_af_plot(parents_v,vcf_sep_comma,nt_to_aa_class, ec19_config['projname'], ec19_lineage, argvs)
        ## reads based analysis
        read_analysis_log = os.path.join(os.path.dirname(argvs.html),'recombinant_reads.log')
        cmd=[os.path.join(bin_dir,'recombinant_read_analysis.py'), '--refacc', 'NC_045512_2','--bam', ec19_bam[0], '-eo', argvs.ec19_projdir, '--vcf', ec19_vcf[0] ,'--igv', argvs.igv]
        if argvs.verbose:
            cmd.append('--verbose')
        process_cmd(cmd,"Running per read analysis", read_analysis_log)
        if 'Omicron' in parents_v and 'Delta' in parents_v:
            deltacron_af_plot_by_sample_id(nt_to_variant, nt_to_aa, delta_uniq_nt,omicron_uniq_nt, ec19_config['projname'], mutations_af, mutations_dp, ec19_lineage, argvs)
    else:
        logging.info(f'No read based recombinant analysis performed since mixed count {mix_count} less than {argvs.minMixed_n} or no parents variants detected.')
         
if __name__ == '__main__':
	main()
