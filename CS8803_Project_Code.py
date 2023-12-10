##### Load Required Libraries #####
import pandas as pd
import numpy as np
import time

##### Helper Functions #####
def residuals(x: int):
    if x == 6:
        return 1
    elif x == -6:
        return -1/3
    elif x == 3:
        return 2/3
    elif x == -3:
        return 0
    else:
        return 1/3

def ttd(start: str, grad: str):
    try:
        years = int(str(start).strip()[:4]) - int(str(grad).strip()[:4])
        remainder = int(str(start).strip()[-2:]) - int(str(grad).strip()[-2:])
        time_to_degree = abs(years) + residuals(remainder)
        return round(time_to_degree, ndigits = 2)
    except:
        return 0.0

def clean_tuple(x: tuple):
    if len(x) > 1:
        return tuple(x)
    else:
        return f"('{x[0]}')"
    
def null_check(x):
    if len(x) == 0:
        return 0
    else:
        return x
    
def letter_grade(x: str):
    x = x.strip('*')
    grade_dictionary = {'A':4, 'B':3, 'C':2, 'D':1, 'F':0}
    try:
        return grade_dictionary[x]
    except:
        return 0.0

def course_avg(x):
    return int(sum(x)/len(x))

def estimate_qlty(x):
    for_gpa = 0
    total_pts = 0
    for grade, count in x['grades'].items():
        value = letter_grade(grade)
        pts = count * value
        total_pts += pts
        if grade in ['A', 'B', 'C', 'D', 'F']:
            for_gpa += count
        else:
            continue
    try:
        return round(total_pts/for_gpa, 2) * x['credits']
    except:
        return 0.0

def gateway_credit(adict: dict):
    for key, value in adict.items():
        if key in ['MATH1113', 'MATH1551', 'MATH1552', 'MATH1553', 'MATH1554', 'MATH2550', 'MATH2551']:
            if value == 1:
                target_courses.add('MATH1113')
                target_courses.add(key)
            else:
                continue
        elif key in ['CS1301', 'CS1371', 'CS1331', 'PHYS2211', 'PHYS2212']:
            if value == 1:
                target_courses.add(key)
            else:
                continue
        else:
            continue

def time_slot(mon: str, tue:str, wed:str, trs:str, fri:str, sat:str, sun:str, start:int, end:int):
    alist = [mon, tue, wed, trs, fri, sat, sun]
    days = [i.upper() for i in alist if len(i) == 1]
    day = "".join(i for i in days)
    return day, start, end

def mod_sql(script, phrase, replacement):
    to_add = clean_tuple(replacement)
    ret_script = script.replace(phrase,str(to_add))
    return ret_script

def check_sets(pool_courses, target_history):
    target_recommendations = pool_courses.difference(target_history)
    unique_targets = target_history.difference(pool_courses)
    return target_recommendations, unique_targets

##### Import and Convert Student Profile Data #####
def build_lake(data, target_id):
    df = pd.read_csv(data)
    df.fillna(0, inplace = True)
    df['TTD'] = df.apply(lambda x: ttd(x['FIRST_TERM'], x['Grad_Term']), axis = 1)
    target_df = df[(df.PIDM == target_id)].copy()
    target_data = target_df.to_dict(orient = 'list')
    gateway_credit(target_data)
    target_df.drop(['FIRST_TERM', 'CURRENT_PROGRAM', 'Changed_Major', 'Grad_Term', 'GPA', 'HRS_EARNED', 'TTD'], axis = 1, inplace = True)
    target_vector = np.array(target_df)
    major = target_data['CURRENT_PROGRAM'][0]
    hours_earned = target_data['HRS_EARNED'][0]
    lake = df[(df['CURRENT_PROGRAM'] == major) & (df['Changed_Major'] == 0) & (df['TTD'] <= 5.1) & ((df['TTD'] > 0) | (df['HRS_EARNED'] >= hours_earned + 30)) & (df['GPA'] >= 3.5)].copy()
    lake.drop(['FIRST_TERM', 'CURRENT_PROGRAM', 'Changed_Major', 'Grad_Term', 'GPA', 'HRS_EARNED', 'TTD'], axis = 1, inplace = True)
    df_array = lake.to_numpy()
    return target_vector, df_array, major, hours_earned

##### Reduce Student Lake to Student Pool #####
def build_pool(target_vector, df_array, threshold = 100):
    pool_set = set()
    loop = 0
    while loop < (target_vector.shape[1] - 1):
        criteria = target_vector.shape[1] - loop
        loop += 1
        target = target_vector[:,1:criteria]
        candidates = df_array[:,1:criteria]
        matches = np.where(np.all(candidates == target,axis = 1))
        pool_set.update(set(matches[0]))
        if len(pool_set) > threshold:
            break
    student_pool = list(df_array[list(pool_set)][:,0])
    if len(student_pool) > 999:
        student_pool = student_pool[:999]
    student_pool.append(target_vector[0][0])
    pool_query = [str(i) for i in student_pool]
    return pool_query, loop

##### Get Course Information for Pool & Target #####
def get_history(pool_query, source_file):
    infile = open(source_file , 'r')
    header = infile.readline()
    data = infile.readlines()
    pool_dict = {}
    target_dict = {}
    course_dict = {}
    course_list = []
    for line in data:
        line_pieces = line.strip().split(',')
        stu_id = int(line_pieces[0])
        term = line_pieces[1]
        crse_subject = line_pieces[2]
        crse_number = line_pieces[3]
        course = str(crse_subject) + str(crse_number)
        crse_crdt = float(null_check(line_pieces[4]))
        alpha_grade = line_pieces[5]
        #title = line_pieces[6]
        num_grade = letter_grade(alpha_grade)
        if str(stu_id) in pool_query:
            if course not in course_dict.keys():
                course_dict[course] = {}
                course_dict[course]['credits'] = []
                course_dict[course]['grades'] = {}
                course_dict[course]['pool'] = set()
                course_dict[course]['pool'].add(stu_id)
                course_dict[course]['credits'].append(crse_crdt)
                course_dict[course]['quality_pts'] = num_grade
                course_dict[course]['grades'][alpha_grade] = 1
                #course_dict[course]['Course Title'] = title
            else:
                course_dict[course]['pool'].add(stu_id)
                course_dict[course]['credits'].append(crse_crdt)
                course_dict[course]['quality_pts'] += num_grade
                if alpha_grade not in course_dict[course]['grades'].keys():
                    course_dict[course]['grades'][alpha_grade] = 1
                else:
                    course_dict[course]['grades'][alpha_grade] += 1
            if stu_id == target_id:
                target_courses.add(course)
                if term not in target_dict.keys():
                    target_dict[term] = {}
                    target_dict[term][course] = alpha_grade
                else:
                    target_dict[term][course] = alpha_grade
            else:
                pool_courses.add(course)
                if stu_id not in pool_dict.keys():
                    pool_dict[stu_id] = {}
                    pool_dict[stu_id][term] = {}
                    pool_dict[stu_id][term][course] = alpha_grade
                else:
                    if term not in pool_dict[stu_id].keys():
                        pool_dict[stu_id][term] = {}
                        pool_dict[stu_id][term][course] = alpha_grade
                    else:
                        pool_dict[stu_id][term][course] = alpha_grade
    for course, detail in course_dict.items():
        newtup = (course, len(course_dict[course]['pool']))
        course_list.append(newtup)
        course_dict[course]['credits'] = course_avg(course_dict[course]['credits'])
        course_dict[course]['estimate_pts'] = estimate_qlty(course_dict[course])
    course_list = sorted(course_list, key=lambda tup: tup[1], reverse = True)
    course_list = [i[0] for i in course_list]
    return target_dict, pool_dict, course_dict, course_list

def term_offerings(source_file):
    infile = open(source_file , 'r')
    header = infile.readline()
    data = infile.readlines()
    term_courses = {}
    for line in data:
        line_pieces = line.strip().split(',')
        crn = int(line_pieces[0])
        subject = str(line_pieces[1])
        crse_num = str(line_pieces[2])
        crdt_hrs = float(null_check(line_pieces[3]))
        seats = int(line_pieces[4])
        start_time = null_check(line_pieces[5])
        end_time = null_check(line_pieces[6])
        mon = str(line_pieces[7])
        tue = str(line_pieces[8])
        wed = str(line_pieces[9])
        trs = str(line_pieces[10])
        fri = str(line_pieces[11])
        sat = str(line_pieces[12])
        sun = str(line_pieces[13])
        seq = str(line_pieces[14])
        pre_sub = str(line_pieces[15])
        seq_num = str(line_pieces[16])
        connector = str(line_pieces[19])
        course_det = subject + crse_num
        times = time_slot(mon, tue, wed, trs, fri, sat, sun, int(start_time), int(end_time))
        if course_det not in term_courses.keys():
            term_courses[course_det] = {}
            term_courses[course_det][crn] = {}
            #term_courses[course_det][crn]['Course'] = subject + crse_num
            term_courses[course_det][crn]['Credits'] = float(crdt_hrs)
            term_courses[course_det][crn]['Seats'] = int(seats)
            term_courses[course_det][crn]['Times'] = set()
            term_courses[course_det][crn]['Times'].add(times)
            term_courses[course_det][crn]['pre_req'] = {}
            term_courses[course_det][crn]['pre_req']['All'] = set()
            term_courses[course_det][crn]['pre_req']['Or'] = set()
            if pre_sub != '':
                if connector == 'A':
                    term_courses[course_det][crn]['pre_req']['All'].add(pre_sub + seq_num)
                elif connector == 'O':
                    term_courses[course_det][crn]['pre_req']['Or'].add(pre_sub + seq_num)
                else:
                    term_courses[course_det][crn]['pre_req']['Or'].add(pre_sub + seq_num) #needs to be adjusted someday
        else:
            if crn not in term_courses[course_det].keys():
                term_courses[course_det][crn] = {}
                term_courses[course_det][crn]['Credits'] = crdt_hrs
                term_courses[course_det][crn]['Seats'] = seats
                term_courses[course_det][crn]['Times'] = set()
                term_courses[course_det][crn]['Times'].add(times)
                term_courses[course_det][crn]['pre_req'] = {}
                term_courses[course_det][crn]['pre_req']['All'] = set()
                term_courses[course_det][crn]['pre_req']['Or'] = set()
                if connector == 'A':
                    term_courses[course_det][crn]['pre_req']['All'].add(pre_sub + seq_num)
                elif connector == 'O':
                    term_courses[course_det][crn]['pre_req']['Or'].add(pre_sub + seq_num)
                else:
                    term_courses[course_det][crn]['pre_req']['Or'].add(pre_sub + seq_num) #need to ammend
            else:
                term_courses[course_det][crn]['Times'].add(times)
                term_courses[course_det][crn]['Seats'] += seats
                if pre_sub != '':
                    if connector == 'A':
                        term_courses[course_det][crn]['pre_req']['All'].add(pre_sub + seq_num)
                    elif connector == 'O':
                        term_courses[course_det][crn]['pre_req']['Or'].add(pre_sub + seq_num)
                    else:
                        term_courses[course_det][crn]['pre_req']['Or'].add(pre_sub + seq_num) #need to ammend
    return term_courses

def prereq_check(target_history, semester_courses, pool_history):
    term_offerings = list(semester_courses.keys())
    qualified_courses = []
    for i in range(len(pool_history)):
        tgt = str(pool_history[i])
        if tgt in term_offerings:
            data = semester_courses[tgt]
            for j in data.keys():
                if j not in qualified_courses:
                    data2 = data[j]
                    mandatory_prereqs = data2['pre_req']['All']
                    opt_prereqs = data2['pre_req']['Or']
                    pre_check = [i for i in mandatory_prereqs if i not in target_history]
                    if len(pre_check) == 0:
                        opt_check = [o for o in opt_prereqs if o in target_history]
                        if len(opt_check) > 0 or len(opt_prereqs) == 0:
                            qualified_courses.append(tgt)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
            else:
                continue
    return qualified_courses

##### Code to conduct test for evaluation #####
dftest = pd.read_csv() #Add file path
dftest = dftest.astype({'PIDM': str})
dftest1 = dftest[(dftest['FIRST_TERM'] == 202208)].copy()
class_2022 = list(dftest1['PIDM'].unique())

def res(x):
    if x == 'Out-of-State Resident*':
        return 'Non_Res'
    else:
        return 'GA_Res'
    
dict_2022 = {}

for i in range(len(class_2022)):
    target_id = int(class_2022[i])
    dict_2022[target_id] = {}
    start_time = time.time()
    threshold = 100
    target_courses = set()
    pool_courses = set()
    target_vector, df_array, major, hrs = build_lake(----, target_id) #Add file path
    pool_query, loop = build_pool(target_vector, df_array)
    target_history, pool_history, course_dictionary, course_list = get_history(pool_query, ----) #Add file path
    target_recommendations, unique_targets = check_sets(pool_courses, target_courses)
    term_courses = term_offerings() #Add file path
    candidate_courses = prereq_check(target_courses, term_courses, course_list)
    total_time = round(time.time() - start_time,2)
    dict_2022[target_id]['Major'] = major
    dict_2022[target_id]['Residency'] = res(target_vector[0][1])
    dict_2022[target_id]['Credits_Earned'] = hrs
    dict_2022[target_id]['Lake_Size'] = len(df_array)
    dict_2022[target_id]['Pool_Size'] = len(pool_query)
    dict_2022[target_id]['Iterations'] = loop
    dict_2022[target_id]['Highest_Rec'] = course_list[0]
    dict_2022[target_id]['Eligible_Courses'] = len(candidate_courses)
    dict_2022[target_id]['Time'] = total_time

df_ret = pd.DataFrame.from_dict(dict_2022, orient = 'index')
#df_ret.to_csv('CS8803_evaluation.csv', encoding = 'utf-8', index = True)


