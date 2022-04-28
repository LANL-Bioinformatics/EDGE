#!/usr/bin/env python3
import os, sys
import glob
import argparse
import re
import plotly.graph_objects as go
from plotly.offline import plot
from plotly.subplots import make_subplots
import json

bin_dir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def setup_argparse():
    parser = argparse.ArgumentParser(description='''Script to do recombinant analysis''')
    parser.add_argument('-m', '--minAF',metavar='[FLOAT]',required=False, type=float, default=0.1, help="plot threshold. minimum average alleic frequency of strains' unqiue mutations [default:0.1]")
    parser.add_argument('--minRecAF',metavar='[FLOAT]',required=False, type=float, default=0.3, help="recobminant threshold. minimum average alleic frequency of strains' unqiue mutations [default:0.3]")
    parser.add_argument('-x', '--mixRatio',metavar='[FLOAT]',required=False, type=float, default=0.5, help="ratio of alleic frequency sites between 0.2 ((minMixAF)) and 0.8 (maxMixAF) to determine the mixed population")
    parser.add_argument('--minMixAF',metavar='[FLOAT]',required=False, type=float, default=0.2, help="minimum alleic frequency for checking mixed mutation [default:0.2]")
    parser.add_argument('--maxMixAF',metavar='[FLOAT]',required=False, type=float, default=0.8, help="maximum alleic frequency for checking mixed mutation [default:0.8]")
    parser.add_argument('-u', '--uniqVar',metavar='[FILE]',required=False, type=str, help="unique variants json file")
    parser.add_argument('-eo', '--ec19_projdir',metavar='[PATH]',required=True, type=str,  help="ec-19 project directory")
    parser.add_argument('--igv', metavar='[PATH]',required=False, type=str,  help="igv.html relative path")
    parser.add_argument('--html' ,metavar='[FILE]',required=False, type=str, help='output plot html')
    argvs = parser.parse_args()

    if not argvs.uniqVar:
        argvs.uniqVar = os.path.join(bin_dir, 'uniq_variants.json')
    if not argvs.html:
        argvs.html = os.path.join(argvs.ec19_projdir, 'ReadsBasedAnalysis','readsMappingToRef','recombinant_analysis_result.html')
    if not argvs.igv:
        argvs.igv = os.path.join('..', '..','IGV','ref_tracks','igv.html')
    return argvs

def parse_variants(vcf,comp,delta_uniq_nt,omicron_uniq_nt,argvs):
    variant_af=dict()
    variant_dp=dict()
    pos_list = [u.split(':')[1] for u in (delta_uniq_nt +  omicron_uniq_nt)]
    comp_pos=dict()
    vcf_comp=dict()
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
            ref_bases = content[3]
            alt_bases = content[4]
            if ',' in alt_bases or (AFreq > argvs.minMixAF and AFreq < argvs.maxMixAF):
                mix_count += 1
            if len(content[3]) - len(content[4]) != 0 and ',' not in alt_bases:
                comp_content = comp_pos[str(int(content[1]) + 1)]
                ref_bases = 'del' if len(content[3]) - len(content[4]) > 0 else 'ins'
                alt_bases = abs(len(content[3]) - len(content[4])) if ref_bases == 'del' else content[4][1:] 
                vcf_comp[ref_bases + ":" + str(int(content[1]) + 1) + ":" + str(alt_bases)] = ':'.join([ comp_content[x].split(' ')[0] for x in range(4,9) ])
            else:
                comp_content = comp_pos[content[1]]
                vcf_comp[ref_bases + ":" + content[1] + ":" + alt_bases] = ':'.join([ comp_content[x].split(' ')[0] for x in range(4,9) ])

    for u in (delta_uniq_nt + omicron_uniq_nt):
        exist = 0
        ref_nt,pos,alt_nt = u.split(':')
   
        for v in vcf_content:
            content = v.strip().split('\t')
            content2 = content[-1].split(':')
            AFreq = content2[-1]
            #total_dp = int(content2[1]) + int(content2[2])
            if ref_nt == 'del' and str(int(pos) - 1) == content[1] and str(len(content[3]) - len(content[4])) == alt_nt:
                #sys.stderr.write ( "\t".join(content) +  str(len(content[3]) - len(content[4])) + "\t" + alt_nt + "\n" )               
                variant_af[u] = AFreq
                exist = 1 
            elif ref_nt == 'ins' and str(int(pos) - 1) == content[1] and alt_nt in content[4]:
                variant_af[u] = AFreq
                exist = 1
            elif pos == content[1] and ref_nt == content[3] and alt_nt == content[4]:
                variant_af[u] = AFreq
                exist = 1
        comp_content = comp_pos[pos]
        all_count = ':'.join([ comp_content[x].split(' ')[0] for x in range(4,9) ])
        variant_dp[u]=str(all_count)
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
                    total_dp = comp_content[3]
                    if ref_base == 'A':
                        count, percentage=comp_content[4].split(' ')
                    if ref_base == 'C':
                        count, percentage=comp_content[5].split(' ')
                    if ref_base == 'G':
                        count, percentage=comp_content[6].split(' ')
                    if ref_base == 'T':
                        count, percentage=comp_content[7].split(' ')
                    percentage = "{:.4f}".format( 1 - (int(count)/int(total_dp))) if int(total_dp) > 0 else '0.000' 

            variant_af[u]=percentage.replace('(','').replace(')','')

   
    return variant_af, variant_dp, vcf_comp, mix_count

def check_variants(variants_af,variants_dp, delta_uniq_nt,omicron_uniq_nt,argvs):
    ## need to set criteria to plot.  default to False.
    check=False
    probably_delta = 0
    probably_omicron = 0
    total_delta_af = 0
    total_omicron_af = 0
    for u in delta_uniq_nt:
        ac,cc,gc,tc,nc = variants_dp[u].split(':')
        total_delta_af = total_delta_af + float(variants_af[u])
        if float(variants_af[u]) >= argvs.maxMixAF:
            probably_delta += 1
        if float(variants_af[u]) <= argvs.minMixAF and (int(ac) + int(cc) + int(gc) + int(tc) + int(nc)) > 0:
            probably_omicron += 1
    for u in omicron_uniq_nt:
        ac,cc,gc,tc,nc = variants_dp[u].split(':')
        total_omicron_af = total_omicron_af + float(variants_af[u])
        if float(variants_af[u]) >= argvs.maxMixAF:
            probably_omicron += 1
        if float(variants_af[u]) <= argvs.minMixAF and (int(ac) + int(cc) + int(gc) + int(tc) + int(nc)) > 0:
            probably_delta += 1
    avg_delta_AF = total_delta_af/len(delta_uniq_nt)
    avg_omicron_AF = total_omicron_af/len(omicron_uniq_nt)
    #if probably_delta > float(minAF) * len(delta_uniq_nt) or probably_omicron > float(minAF) * len(omicron_uniq_nt):
    if avg_delta_AF > argvs.minAF or avg_omicron_AF > argvs.minAF:
        check=True
    sys.stderr.write(f'Average Delta Unique Variants AF: {avg_delta_AF}\n')
    sys.stderr.write(f'Average Omicron Unique Variants AF: {avg_omicron_AF}\n')
    sys.stderr.write(f'Probably Delta Unique Variants Count (>={argvs.maxMixAF}): {probably_delta}\n')
    sys.stderr.write(f'Probably Omicron Unique Variants Count (>={argvs.maxMixAF}): {probably_omicron}\n')
    if avg_delta_AF > argvs.minRecAF and avg_omicron_AF > argvs.minRecAF and int(len(delta_uniq_nt) * argvs.minRecAF) < probably_delta  and int(len(omicron_uniq_nt) * argvs.minRecAF) < probably_omicron :
        sys.stderr.write(f'Probable Delta and Omicron recombinant. [minimum average delta and omicron unqiue variants AF: {argvs.minRecAF}]\n')
    return check

def load_uniq_var(file):
    with open( file, 'r') as f:
        unique_v_data = json.load(f)
    delta_uniq_nt=list(unique_v_data['delta'].keys())
    delta_uniq_aa= list(unique_v_data['delta'].values())
    omicron_uniq_nt=list(unique_v_data['omicron'].keys())
    omicron_uniq_aa=list(unique_v_data['omicron'].values())
    nt_to_aa = {**unique_v_data['delta'],**unique_v_data['omicron']}
    new_d = { i: int(i.split(':')[1]) for i in nt_to_aa }
    variants_nt = list(dict(sorted(new_d.items(), key=lambda item: item[1])).keys())
    variants_aa = [ nt_to_aa[x] for x in variants_nt ]
    return(delta_uniq_nt,omicron_uniq_nt,variants_nt,variants_aa)

def genome_af_plot(variants,variants_aa,delta_uniq_nt,omicron_uniq_nt,variants_af,barplot_lists, lineage , sample, url , output):
    fig = make_subplots(rows=2,cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[200,500])
   
    fig.add_trace(go.Scatter(
        x= variants,
        y=[ float(variants_af[i]) if i in delta_uniq_nt else None for i in variants ],
        mode="markers",
        marker=dict(color='red',size=8),
        name="Delta unique variants",
        hovertemplate = 'AF: %{y:.2f}<br>'+'<b>%{text}</b>',
        text=[ f'{i}' for i in variants_aa ],
        customdata=[ url + '?locus=NC_045512_2:' + str(int(i.split(':')[1]) - 100) + '-' +  str(int(i.split(':')[1]) + 100) if i in delta_uniq_nt else None for i in variants ],
        showlegend=True,
        legendgroup=1
    ),row=2,col=1)
    fig.add_trace(go.Scatter(
        x=variants,
        y=[ float(variants_af[i]) if i in omicron_uniq_nt else None for i in variants ],
        mode="markers",
        marker=dict(color='blue',size=8),
        name="Omicron unique variants",
        hovertemplate = 'AF: %{y:.2f}<br>'+'<b>%{text}</b>',
        text=[ f'{i}' for i in variants_aa ],
        customdata=[ url + '?locus=NC_045512_2:' + str(int(i.split(':')[1]) - 100) + '-' +  str(int(i.split(':')[1]) + 100) if i in omicron_uniq_nt else None for i in variants ],
        showlegend=True,
        legendgroup=1
    ),row=2,col=1)
    fig.add_trace(go.Bar(
        x=variants,
        y=barplot_lists[0],
        marker_color='red',
        name="Delta Depth Coverage",
        showlegend=False,
        offsetgroup=0,
    ),row=1,col=1)
    fig.add_trace(go.Bar(
        x=variants,
        y=barplot_lists[1],
        marker_color='blue',
        name="Omicron Depth Coverage",
        showlegend=False,
        offsetgroup=0,
        base=barplot_lists[0]
    ),row=1,col=1)
    fig.add_trace(go.Bar(
        x=variants,
        y=barplot_lists[2],
        marker_color='grey',
        name="Other Depth Coverage",
        showlegend=False,
        offsetgroup=0,
        base=[val1+val2 for val1, val2 in zip(barplot_lists[0],barplot_lists[1])],
    ),row=1,col=1)
    fig.add_vline(x=10.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=14.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=40.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=41.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=42.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=45.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=47.5, line_width=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=48.5, line_width=1, line_dash="dash", line_color="grey")
            
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

def bar_plot_list3(variants, delta_uniq_nt,omicron_uniq_nt,variants_dp):
    count_dict=dict()
    delta_count_list=list()
    omicron_count_list=list()
    other_count_list=list()
    for i in variants:
        ref_nt,pos,alt_nt = i.split(':')
        ac,cc,gc,tc,nc = variants_dp[i].split(':')
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

def genome_af_plot_by_sample_id(variants_nt,variants_aa,delta_uniq_nt,omicron_uniq_nt, sample, variants_af, variants_dp, lineage, igv_url , out_html):

    barplot_lists = bar_plot_list3(variants_nt,delta_uniq_nt,omicron_uniq_nt, variants_dp)
    genome_af_plot(variants_nt,variants_aa,delta_uniq_nt,omicron_uniq_nt,variants_af, barplot_lists, lineage, sample, igv_url,out_html)

def comp_stack_bar_plot(vcf_comp, mix_count, delta_uniq_nt,omicron_uniq_nt,  argvs):
    barplot_html_output = os.path.join(os.path.dirname(argvs.html),'variants_bar_plot.html')
    count_dict=dict()
    vcf_list = list(vcf_comp.keys())
    aclist, cclist, gclist, tclist, nclist = list(), list(),list(),list(),list()
    aplist, cplist, gplist, tplist, nplist = list(), list(),list(),list(),list()
    mutations_count = len(vcf_comp)
    mix_ratio = mix_count/mutations_count * 100 if mutations_count > 0 else 0
    for i in vcf_comp:
        ac,cc,gc,tc,nc = vcf_comp[i].split(':')
        total = int(ac) + int(cc) + int(gc) + int(tc) + int(nc)
        aclist.append(int(ac))
        aplist.append(int(ac)/total)
        cclist.append(int(cc))
        cplist.append(int(cc)/total)
        gclist.append(int(gc))
        gplist.append(int(gc)/total)
        tclist.append(int(tc))
        tplist.append(int(tc)/total)
        nclist.append(int(nc))
        nplist.append(int(nc)/total)
    fig = make_subplots(rows=3,cols=1, shared_xaxes=True, vertical_spacing=0.01, row_heights=[200,500,15])
    fig.add_trace(go.Bar(name='A', x=vcf_list, y=aclist,marker_color='green',opacity=.7),row=1,col=1)
    fig.add_trace(go.Bar(name='C', x=vcf_list, y=cclist,marker_color='blue',opacity=.7),row=1,col=1)
    fig.add_trace(go.Bar(name='G', x=vcf_list, y=gclist,marker_color='orange',opacity=.7),row=1,col=1)
    fig.add_trace(go.Bar(name='T', x=vcf_list, y=tclist,marker_color='red',opacity=.7),row=1,col=1)
    fig.add_trace(go.Bar(name='INDELs/N', x=vcf_list, y=nclist,marker_color='grey',opacity=.7),row=1,col=1)

    fig.add_trace(go.Bar(name='A', x=vcf_list, y=aplist,marker_color='green',showlegend=False,opacity=.7),row=2,col=1)
    fig.add_trace(go.Bar(name='C', x=vcf_list, y=cplist,marker_color='blue', showlegend=False,opacity=.7),row=2,col=1)
    fig.add_trace(go.Bar(name='G', x=vcf_list, y=gplist,marker_color='orange', showlegend=False,opacity=.7),row=2,col=1)
    fig.add_trace(go.Bar(name='T', x=vcf_list, y=tplist,marker_color='red', showlegend=False,opacity=.7),row=2,col=1)
    fig.add_trace(go.Bar(name='INDELs/N', x=vcf_list, y=nplist,marker_color='grey', showlegend=False,opacity=.7),row=2,col=1)

    omicron_list = [ 1 if i in omicron_uniq_nt else None for i in vcf_comp ]
    delta_list = [ 1 if i in delta_uniq_nt else None for i in vcf_comp ]
    fig.add_trace(go.Bar(name='Omicron', x=vcf_list, y=omicron_list,marker_color='blue', hovertemplate="Omicron", showlegend=False,opacity=.3),row=3,col=1)
    fig.add_trace(go.Bar(name='Delta', x=vcf_list, y=delta_list,marker_color='red', hovertemplate="Delta", showlegend=False,opacity=.3),row=3,col=1)

    fig.update_layout(barmode='stack', title_text='Nucleotide Composition of Mutations')
    fig.update_xaxes(tickfont=dict(size=8),tickangle=-60)
    fig.update_yaxes(title="Count",row=1,col=1)
    fig.update_yaxes(title="Percentage",row=2,col=1)
    fig.update_yaxes(visible=False,row=3,col=1)
    if mix_count > 0:
        fig.add_annotation(text=f"{mix_ratio:.2f}% ({mix_count}/{mutations_count}) mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF}.",
                  xref="x domain", yref="y domain", showarrow=False, font=dict(size=10),
                  x=0, y=1.12,row=1,col=1)
    fig.write_html(barplot_html_output)
    fig.write_image(barplot_html_output+'.png')      

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

def main():
    argvs = setup_argparse()
    (delta_uniq_nt,omicron_uniq_nt,variants_nt,variants_aa) = load_uniq_var(argvs.uniqVar)

    ec19_consensus = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.fasta",recursive = False)
    ec19_log = glob.glob(argvs.ec19_projdir+f"/process.log",recursive = False)
    ec19_config_file = glob.glob(argvs.ec19_projdir+f"/config.txt",recursive = False)
    ec19_lineage_file = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus_lineage.txt",recursive = False)
    ec19_vcf = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.vcf",recursive = False)
    ec19_snp_result = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.SNPs_report.txt",recursive = False)
    ec19_indel_result = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.Indles_report.txt",recursive = False)
    ec19_compositionlog = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2_consensus.compositionlog",recursive = False)
    ec19_aln_stats = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/readsMappingToRef/NC_045512.2.alnstats.txt",recursive = False)
    ec19_fq_count = glob.glob(argvs.ec19_projdir+f"/QcReads/fastqCount.txt",recursive = False)
    ec19_lineage_abund = glob.glob(argvs.ec19_projdir+f"/ReadsBasedAnalysis/LineageAbundance/predicitions.tsv",recursive = False)
    igv_relative_url = argvs.igv
    # cdc_data/aegis/2021-12-20/313355412-1639983320/ASC210517005-B2968347/HTML_Report/igv.html
    
    if not ec19_vcf:
        sys.stderr.write("Cannot find consensus.vcf file\n")
        exit
    if not ec19_compositionlog:
        sys.stderr.write("Cannot find consensus.compositionlog\n")
        exit

    ec19_lineage = parse_lineage(ec19_lineage_file[0]) if ec19_lineage_file else 'Unknown'
    ec19_config = parse_ec19_config(ec19_config_file[0])
    variants_af, variants_dp, vcf_comp, mix_count = parse_variants(ec19_vcf[0],ec19_compositionlog[0],delta_uniq_nt,omicron_uniq_nt,argvs)
    comp_stack_bar_plot(vcf_comp,mix_count,delta_uniq_nt,omicron_uniq_nt, argvs)
    plot_bool = check_variants(variants_af,variants_dp, delta_uniq_nt,omicron_uniq_nt,argvs)
    if plot_bool:
        genome_af_plot_by_sample_id(variants_nt,variants_aa,delta_uniq_nt,omicron_uniq_nt, ec19_config['projname'], variants_af, variants_dp, ec19_lineage, igv_relative_url, argvs.html)
        
    mix_ratio = mix_count/len(vcf_comp) * 100 if len(vcf_comp) > 0 else 0
    mutations_count = len(vcf_comp)
    sys.stderr.write(f"({mix_count}/{mutations_count}) mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF}. \n")
    if len(vcf_comp) > 0 and mix_count/len(vcf_comp) > argvs.mixRatio:
        sys.stderr.write(f"Probable Mixed Infection. {mix_ratio:.2f}% ({mix_count}/{mutations_count}) mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF}. \n")
    
if __name__ == '__main__':
	main()