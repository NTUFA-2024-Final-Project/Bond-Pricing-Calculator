import React, { useState } from 'react';
import PropTypes from 'prop-types';
import '../App/App.css';
import axios from 'axios';
import Select from 'react-select';

const initialValues = {
	/* -- basic information -- */
    k: '10.00',
	maturity_date: '2030-05-31',
	delta: '0.98',
	pi: '0.50',

	/* -- Callable Bond Information -- */
	coupon_rate: '0.10',
	face_value: '50.00',
	day_count_convention:'0', 

}

const BondForm = ({ change }) => {
	const [state, setState] = useState(initialValues);

    // Options for the day count convention dropdown
    const dayCountOptions = [
        { value: '0', label: 'Thirty/360' },
        { value: '1', label: 'Thirty/365' },
        { value: '2', label: 'Actual/360' },
        { value: '3', label: 'Actual/365' },
        { value: '4', label: 'Actual/Actual' }
    ];

	const handleChange = e => {
		let { value, name, type } = e.target;

		if (type === 'file') {
			setState(prevState => ({
				...prevState,
				[name]: e.target.files[0] // Capture the file itself
			}));
		} else {
			// Handling other inputs
			setState(prevState => ({
				...prevState,
				[name]: value
			}));
		}	
	};

	const handleSelectChange = (selectedOption) => {
        setState(prev => ({ ...prev, day_count_convention: selectedOption }));
    };

	const handleSubmit = () => {
		change(state);
		setState(initialValues);
	};

	// const handleSubmit = async () => {
	// 	const formData = new FormData();
		
	// 	// Append file to FormData if present
	// 	if (state.upload_file) {
	// 		formData.append('upload_file', state.upload_file);
	// 	}
	
	// 	// Append other form data as JSON
	// 	const json = JSON.stringify({
	// 		k: state.k,
	// 		maturity_date: state.maturity_date,
	// 		delta: state.delta,
	// 		pi: state.pi,
	// 		coupon_rate: state.coupon_rate,
	// 		face_value: state.face_value,
	// 		day_count_convention: state.day_count_convention
	// 	});
	// 	formData.append('json', new Blob([json], {
	// 		type: 'application/json'
	// 	}));
	
	// 	try {
	// 		const response = await axios.post('http://localhost:3001/calculate', formData, {
	// 			headers: {
	// 				'Content-Type': 'multipart/form-data'
	// 			}
	// 		});
	// 		console.log('Submission successful', response.data);
	// 	} catch (error) {
	// 		console.error('Error submitting form:', error);
	// 	}
	// };
	
	return (
		<>
			<div className="row">
				<div className='row center'>
					<h4 className='white-text'> Basic Information </h4>
				</div>

				<div className="center col s3">
					<label htmlFor="k"> Strike Price </label>
					<input
						id="k"
						name="k"
						type="number"
						min="0.00"
						max="200.00"
						step="0.01"
						placeholder="100.00"
						value={state.k}
						onChange={handleChange}
					/>
				</div>

				<div className="center col s3">
					<label htmlFor="maturity_date"> Maturity Date </label>
					<input
						id="maturity_date"
						name="maturity_date"
						type="date"
						value={state.maturity_date}
						onChange={handleChange}
					/>
				</div>

				<div className="center col s3">
					<label htmlFor="delta"> Delta </label>
					<input
						id="delta"
						name="delta"
						type="number"
						min="0.00"
						max="100.00"
						step="0.01"
						placeholder="0.98"
						value={state.delta}
						onChange={handleChange}
						readOnly
					/>
				</div>

				<div className="center col s3">
					<label htmlFor="pi">Pi </label>
					<input
						id="pi"
						name="pi"
						type="number"
						min="0"
						max="100.00"
						placeholder="0.50"
						value={state.pi}
						onChange={handleChange}
						readOnly
					/>
				</div>
				
				<br/><br/><br/><br/>
				<div className='row center'>
					<h4 className='white-text'> Callable Bond Information </h4>
				</div>
				
				<div className="center col s3">
					<label htmlFor="coupon_rate">Coupon Rate (as Decimal) </label>
					<input
						id="coupon_rate"
						name="coupon_rate"
						type="number"
						min="0"
						max="1.00"
						step="0.01"
						placeholder="0.10"
						value={state.coupon_rate}
						onChange={handleChange}
					/>
				</div>
				
				<div className="center col s3">
					<label htmlFor="face_value">Face Value </label>
					<input
						id="face_value"
						name="face_value"
						type="number"
						min="0"
						max="100.00"
						placeholder="50.00"
						value={state.face_value}
						onChange={handleChange}
					/>
				</div>

				<div className="center col s3">
					<label htmlFor="day_count_convention">Day Count Convention</label>
					<Select
						id="day_count_convention"
						name="day_count_convention"
						value={state.day_count_convention}
						onChange={handleSelectChange}
						options={dayCountOptions}
					/>
				</div>
			</div>
			
			<div className="center">
				<button
					id="bond-btn"
					className="calculate-btn"
					type="button"
					disabled={state.k  === '' || 
							  state.maturity_date === '' || 
							  state.delta === '' || 
							  state.pi === '' || 
							  state.coupon_rate === '' || 
							  state.face_value === '' || 
							  state.day_count_convention === ''}
					onClick={handleSubmit}
				>
					Calculate Bond Price
				</button>
			</div>
		</>
	);
};

BondForm.propTypes = {
	change: PropTypes.func.isRequired
};

export default BondForm;
