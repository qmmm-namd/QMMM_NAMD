#!/usr/bin/env python3 
import numpy as np 
from tools_qmmm import unit
import copy 

class array:
    def __init__(self, inp, group, LA_list, CS_list, \
                qm_array, mm_array, qm_mm_region, b_info_list, \
                    QM_LA_region, ref_dic, grad=True):
        self.inp = inp 
        self.grad = grad 

        self.qm_array = qm_array
        self.mm_array = mm_array
        
        self.QM_LA_region = np.array(QM_LA_region)
        self.ref_dic = ref_dic 
        
        self.group = group
        # [q1, m1]
        self.LA_list = LA_list
        # [m1, m2]
        self.CS_list = CS_list

        self.qm_mm_region = np.array(qm_mm_region)
        
        self.b_info_list = b_info_list
        assert len(self.b_info_list) == 3

        return 


    def adjust_array(self):
        # print(*self.ref_dic.items(), sep='\n')
        # reshape the complete size with canceling terms of link atom and redistributed point charge
        self.reshape_region()

        if self.grad == 1:
            # evaluate bonded gradient with LA, and cancel the corresponding term 
            self.get_la_grad()

        # remove link atoms and correct the gradient 
        self.remove_LA()
        # remove the redistributed point charge and correct the gradient  
        self.remove_rd_charge()
            
        # update array term to overall region 
        self.update_array()
        
        return self.qm_array, self.mm_array, self.LA_array
    
    
    def reshape_region(self):
        if len(self.LA_list) > 0:
            self.LA_array = self.qm_array[-len(self.LA_list):]
            del self.qm_array[-len(self.LA_list):]
            
        if len(self.CS_list) > 0:
            self.CS_array = self.mm_array[-len(self.CS_list):]
            del self.mm_array[-len(self.CS_list):]
        
        assert type(self.qm_array) == list and type(self.mm_array) == list
        self.qm_mm_array = [x[0] for x in sorted(zip(self.qm_array + self.mm_array, self.group.QM + self.group.MM), \
            key = lambda x:x[1])]
        self.qm_mm_array = np.array(self.qm_mm_array)
        
        return 
    

    def update_array(self):
        self.qm_array = [self.qm_mm_array[x-1] for x in self.group.QM]
        self.mm_array = [self.qm_mm_array[x-1] for x in self.group.MM]
        return 


    def get_la_grad(self):
        la_grad = np.zeros(self.QM_LA_region.shape)
        la_grad += self.bond_grad()
        la_grad += self.angle_grad()
        la_grad += self.dihedral_grad()
        
        new_grad = np.zeros(la_grad.shape)
        if len(self.ref_dic) > 0:
            for i in self.ref_dic:
                qm_mm_index = int(i) - 1
                qm_la_index = int(self.ref_dic[i]) - 1
                new_grad[qm_mm_index] = la_grad[qm_la_index]
        else:
            new_grad = copy.deepcopy(la_grad)
                
        for index, term in enumerate(self.group.QM):
            self.qm_mm_array[term-1] -= new_grad[index]
        
        la_array = new_grad[-len(self.LA_list):]
        self.LA_array = np.array(self.LA_array)
        self.LA_array -= la_array
        self.LA_array = self.LA_array.tolist()
        
        return self.LA_array
    
    

    def remove_LA(self):
        # distance list for every atom pair in LA_list, a fixed distance 
        dis_list = [self.inp.LA_dis * unit.ang_2_bohr] * len(self.LA_list)

        self.correct_array(array_list=self.LA_array, index_list=self.LA_list, dis_list=dis_list)
        return 


    def remove_rd_charge(self):
        if len(self.CS_list) == 0:
            return 
        
        # distance list for every atom pair in CS_list, half distance between atom pair
        dis_list = []
        for m1, m2 in self.CS_list:
            m1_coor = self.qm_mm_region[m1-1]
            m2_coor = self.qm_mm_region[m2-1]
            dis = ((m1_coor - m2_coor) ** 2).sum() ** 0.5 / 2 
            dis_list.append(dis)
            
        self.correct_array(array_list=self.CS_array, index_list=self.CS_list, dis_list=dis_list, CS=True)
        
        return 

    
    def correct_array(self, array_list, index_list, dis_list, CS=False):
        for index, i in enumerate(index_list):
            dis = dis_list[index]
            # q1 and m1 are not especially defined as quantum atom and mechanical atom, but two different atoms 
            q1, m1 = i 
            # assert q1 in self.group.QM and m1 in self.group.MM
            q1_coor = self.qm_mm_region[q1-1]
            m1_coor = self.qm_mm_region[m1-1]
            # dis_q1_m1 = bondedtool.cal_distance(coor_list=[q1_coor.tolist(), m1_coor.tolist()])
            dis_q1_m1 = ((q1_coor - m1_coor) ** 2).sum() ** 0.5

            la_array = np.array(array_list[index])

            if CS:
                self.qm_mm_array[q1-1] += 0.5 * la_array
                self.qm_mm_array[m1-1] += 0.5 * la_array

                
            else:
                # chain rule 
                g = dis / dis_q1_m1
                # normal vector of q1-m1
                n = (m1_coor - q1_coor) / np.linalg.norm((m1_coor - q1_coor))
                self.qm_mm_array[m1-1] += g * la_array - g * np.dot(la_array, n) * n
                self.qm_mm_array[q1-1] += (1 - g) * la_array + g * np.dot(la_array, n) * n
                
                # # Walker, R.C., Crowley, M.F. & Case, D.A. J. Comput. Chem. 29, 1019–1031 (2008). 
                # d_m1_x = dis / dis_q1_m1 * np.ones(q1_coor.shape) - dis / dis_q1_m1 ** 3 * (m1_coor - q1_coor) * (m1_coor - q1_coor)[0]
                # d_m1_y = dis / dis_q1_m1 * np.ones(q1_coor.shape) - dis / dis_q1_m1 ** 3 * (m1_coor - q1_coor) * (m1_coor - q1_coor)[1]
                # d_m1_z = dis / dis_q1_m1 * np.ones(q1_coor.shape) - dis / dis_q1_m1 ** 3 * (m1_coor - q1_coor) * (m1_coor - q1_coor)[2]
                
                # d_q1_x = (1 - dis / dis_q1_m1) * np.ones(q1_coor.shape) + dis / dis_q1_m1 ** 3 * (m1_coor - q1_coor) * (m1_coor - q1_coor)[0]
                # d_q1_y = (1 - dis / dis_q1_m1) * np.ones(q1_coor.shape) + dis / dis_q1_m1 ** 3 * (m1_coor - q1_coor) * (m1_coor - q1_coor)[1]
                # d_q1_z = (1 - dis / dis_q1_m1) * np.ones(q1_coor.shape) + dis / dis_q1_m1 ** 3 * (m1_coor - q1_coor) * (m1_coor - q1_coor)[2]
                
                # self.qm_mm_array[m1-1][0] += np.dot(d_m1_x, la_array)
                # self.qm_mm_array[m1-1][1] += np.dot(d_m1_y, la_array)
                # self.qm_mm_array[m1-1][2] += np.dot(d_m1_z, la_array)
                
                # self.qm_mm_array[q1-1][0] += np.dot(d_q1_x, la_array)
                # self.qm_mm_array[q1-1][1] += np.dot(d_q1_y, la_array)
                # self.qm_mm_array[q1-1][2] += np.dot(d_q1_z, la_array)

        return 




    def dihedral_grad(self):
        # change unit of coordinates to angstrom with respect to constants
        # qm_mm_region = self.qm_mm_region / unit.ang_2_bohr
        qm_mm_region = self.QM_LA_region / unit.ang_2_bohr
        f = np.zeros(qm_mm_region.shape)
        for i in self.b_info_list[2]:
            mi, mj, mk, ml = i[-1]
            ci = qm_mm_region[mi-1]
            cj = qm_mm_region[mj-1]
            ck = qm_mm_region[mk-1]
            cl = qm_mm_region[ml-1]
            force_constant = i[0]
            phase = i[1]
            # phase shift in radians
            pn = i[2]
            # curr_value = self.cal_dihedral([ci, cj, ck, cl])
            
            cij = cj - ci 
            cjk = ck - cj 
            ckl = cl - ck 
            
            n1 = np.cross(cij, cjk)
            n2 = np.cross(cjk, ckl)
            n1n = np.linalg.norm(n1)
            n2n = np.linalg.norm(n2)
            
            dot = np.dot(n1, n2)
            d = dot / (n1n * n2n)
            
            curr_value = np.arccos(d)
            
            dn1xi = (n1[1] * cjk[2] - n1[2] * cjk[1]) / n1n 
            dn1xk = (n1[1] * cij[2] - n1[2] * cij[1]) / n1n 
            dn1xj = -1 * (dn1xi + dn1xk)
            dn1xl = 0
            
            dn1yi = (n1[2] * cjk[0] - n1[0] * cjk[2]) / n1n
            dn1yk = (n1[2] * cij[0] - n1[0] * cij[2]) / n1n 
            dn1yj = -1 * (dn1yi + dn1yk)
            dn1yl = 0
            
            dn1zi = (n1[0] * cjk[1] - n1[1] * cjk[0]) / n1n 
            dn1zk = (n1[0] * cij[1] - n1[1] * cij[0]) / n1n 
            dn1zj = -1 * (dn1zi + dn1zk)
            dn1zl = 0
            

            dn2xi = 0
            dn2xj = (n2[1] * ckl[2] - n2[2] * ckl[1]) / n2n 
            dn2xl = (n2[1] * cjk[2] - n2[2] * cjk[1]) / n2n 
            dn2xk = -1 * (dn2xj + dn2xl)
            
            dn2yi = 0
            dn2yj = (n2[2] * ckl[0] - n2[0] * ckl[2]) / n2n 
            dn2yl = (n2[2] * cjk[0] - n2[0] * cjk[2]) / n2n 
            dn2yk = -1 * (dn2yj + dn2yl)
            
            dn2zi = 0
            dn2zj = (n2[0] * ckl[1] - n2[1] * ckl[0]) / n2n 
            dn2zl = (n2[0] * cjk[1] - n2[1] * cjk[0]) / n2n 
            dn2zk = -1 * (dn2zj + dn2zl)
            
            
            ddotxi = n2[1] * cjk[2] - n2[2] * cjk[1]
            ddotyi = n2[2] * cjk[0] - n2[0] * cjk[2] 
            ddotzi = n2[0] * cjk[1] - n2[1] * cjk[0]
            
            ddotxl = n1[1] * cjk[2] - n1[2] * cjk[1] 
            ddotyl = n1[2] * cjk[0] - n1[0] * cjk[2] 
            ddotzl = n1[0] * cjk[1] - n1[1] * cjk[0] 
            
            ddotxj = n2[2] * (cjk[1] + cij[1]) - n2[1] * (cij[2] + cjk[2]) + n1[1] * ckl[2] - n1[2] * ckl[1] 
            ddotyj = n2[0] * (cjk[2] + cij[2]) - n2[2] * (cij[0] + cjk[0]) + n1[2] * ckl[0] - n1[0] * ckl[2] 
            ddotzj = n2[1] * (cjk[0] + cij[0]) - n2[0] * (cij[1] + cjk[1]) + n1[0] * ckl[1] - n1[1] * ckl[0] 
            
            ddotxk = n1[2] * (ckl[1] + cjk[1]) - n1[1] * (cjk[2] + ckl[2]) + n2[1] * cij[2] - n2[2] * cij[1] 
            ddotyk = n1[0] * (ckl[2] + cjk[2]) - n1[2] * (cjk[0] + ckl[0]) + n2[2] * cij[0] - n2[0] * cij[2] 
            ddotzk = n1[1] * (ckl[0] + cjk[0]) - n1[0] * (cjk[1] + ckl[1]) + n2[0] * cij[1] - n2[1] * cij[0] 
            
        
            dn1i = np.array([dn1xi, dn1yi, dn1zi])
            dn1j = np.array([dn1xj, dn1yj, dn1zj])
            dn1k = np.array([dn1xk, dn1yk, dn1zk])
            dn1l = np.array([dn1xl, dn1yl, dn1zl])
            
            dn2i = np.array([dn2xi, dn2yi, dn2zi])
            dn2j = np.array([dn2xj, dn2yj, dn2zj])
            dn2k = np.array([dn2xk, dn2yk, dn2zk])
            dn2l = np.array([dn2xl, dn2yl, dn2zl])
            
            ddoti = np.array([ddotxi, ddotyi, ddotzi])
            ddotj = np.array([ddotxj, ddotyj, ddotzj])
            ddotk = np.array([ddotxk, ddotyk, ddotzk])
            ddotl = np.array([ddotxl, ddotyl, ddotzl])
            
            ddi = (ddoti * n1n * n2n - (dn1i * n2n + n1n * dn2i) * dot) / ((n1n * n2n) ** 2 )
            ddj = (ddotj * n1n * n2n - (dn1j * n2n + n1n * dn2j) * dot) / ((n1n * n2n) ** 2 )
            ddk = (ddotk * n1n * n2n - (dn1k * n2n + n1n * dn2k) * dot) / ((n1n * n2n) ** 2 )
            ddl = (ddotl * n1n * n2n - (dn1l * n2n + n1n * dn2l) * dot) / ((n1n * n2n) ** 2 )
          
            dphii = -1 / ((1 - d ** 2) ** 0.5) * ddi
            dphij = -1 / ((1 - d ** 2) ** 0.5) * ddj
            dphik = -1 / ((1 - d ** 2) ** 0.5) * ddk
            dphil = -1 / ((1 - d ** 2) ** 0.5) * ddl
            
            # # only correct forces in QM region 
            # if mi in self.group.QM:
            #     f[mi-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphii  
            # if mj in self.group.QM:
            #     f[mj-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphij 
            # if mk in self.group.QM:
            #     f[mk-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphik 
            # if ml in self.group.QM:
            #     f[ml-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphil 
            
            # # only correct forces in MM region 
            # if mi in self.group.MM:
            #     f[mi-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphii  
            # if mj in self.group.MM:
            #     f[mj-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphij 
            # if mk in self.group.MM:
            #     f[mk-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphik 
            # if ml in self.group.MM:
            #     f[ml-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphil 

            f[mi-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphii  
            f[mj-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphij 
            f[mk-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphik 
            f[ml-1] += force_constant * pn * np.sin(pn * curr_value - phase) * dphil 

        # force to gradient, and kcal/angstrom to au
        g = -1 * f * unit.kcal_2_au / unit.ang_2_bohr
        # self.update_la_array(g)

        return g
       

    def angle_grad(self):
        # change unit of coordinates to angstrom with respect to constants
        # qm_mm_region = self.qm_mm_region / unit.ang_2_bohr
        qm_mm_region = self.QM_LA_region / unit.ang_2_bohr
        f = np.zeros(qm_mm_region.shape)
        for i in self.b_info_list[1]:
            mi, mj, mk = i[-1]
            # atom i, j, k. j is the middle one
            ci = qm_mm_region[mi-1]
            cj = qm_mm_region[mj-1]
            ck = qm_mm_region[mk-1]
            cij = cj - ci 
            ckj = cj - ck
            
            rij = ((ci - cj) ** 2).sum() ** 0.5 
            rkj = ((cj - ck) ** 2).sum() ** 0.5 
            
            ddi = (-rij * ckj + cij / rij * np.dot(cij, ckj)) / (rij ** 2 * rkj)
            ddk = (-cij * rkj + ckj / rkj * np.dot(cij, ckj)) / (rij * rkj ** 2)
            
            delta = np.dot(cij, ckj) / (rij * rkj)
            phi = np.degrees(np.arccos(delta))
            
            
            di = - 1 / ((1 - delta ** 2) ** 0.5) * ddi
            dk = - 1 / ((1 - delta ** 2) ** 0.5) * ddk
            
            k = i[0]
            ref_ = i[1]
            
            fi = - 2 * k * np.deg2rad(phi - ref_) * di
            fk = - 2 * k * np.deg2rad(phi - ref_) * dk
            
            fj = -(fi + fk)
            
            # # only correct forces in QM region 
            # if mi in self.group.QM:
            #     f[mi-1] += fi
            # if mj in self.group.QM:
            #     f[mj-1] += fj
            # if mk in self.group.QM:
            #     f[mk-1] += fk
            
            # # only correct forces in MM region 
            # if mi in self.group.MM:
            #     f[mi-1] += fi
            # if mj in self.group.MM:
            #     f[mj-1] += fj
            # if mk in self.group.MM:
            #     f[mk-1] += fk


            f[mi-1] += fi
            f[mj-1] += fj
            f[mk-1] += fk
            
        # force to gradient and kcal/angstrom to au 
        g = -1 * f * unit.kcal_2_au / unit.ang_2_bohr
        # self.update_la_array(g)
        # self.qm_mm_array -= g
        
        return g
    
    
    def bond_grad(self):
        # change unit of coordinates to angstrom with respect to constants
        # qm_mm_region = self.qm_mm_region / unit.ang_2_bohr
        qm_mm_region = self.QM_LA_region / unit.ang_2_bohr
        f = np.zeros(qm_mm_region.shape)
        for i in self.b_info_list[0]:
            m1, m2 = i[-1]
            c1, c2 = qm_mm_region[m1-1], qm_mm_region[m2-1]
            rk = i[0]
            req = i[1]
            curr_value = ((c1 - c2) ** 2).sum() ** 0.5 
            f_tmp = 2 * rk * (curr_value - req) * (c1 - c2) / curr_value 

            # only correct forces in QM region 
            # if m1 in self.group.QM:
            #     f[m1-1] -= f_tmp
            # if m2 in self.group.QM:
            #     f[m2-1] += f_tmp
            
            # # only correct forces in <M region 
            # if m1 in self.group.MM:
            #     f[m1-1] -= f_tmp
            # if m2 in self.group.MM:
            #     f[m2-1] += f_tmp


            f[m1-1] -= f_tmp
            f[m2-1] += f_tmp
             
        # force to gradient and kcal/angstrom to au    
        g = -1 * f * unit.kcal_2_au / unit.ang_2_bohr
        # self.update_la_array(g)
        # self.qm_mm_array -= g
        
        return g
    
